function [reconstructed, residual] = ft_singletrialanalysis(cfg, data)

% FT_SINGLETRIALANALYSIS computes a single-trial estimate of the event-
% related activity
%
% Use as
%   [reconstructed, residual] = ft_aseo(cfg, data)
% where data is single-channel raw data as obtained by FT_PREPROCESSING
% and cfg is a configuration structure according to
%
%  cfg.method                       = different methods of estimating event-related activity
%                                        'aseo', analysis of single-trial ERP and ongoing
%                                        activity (according to Xu et al, 2009)
%  cfg.channel                      = Nx1 cell-array with selection of channels (default = 'all'),
%                                       see FT_CHANNELSELECTION for details
%  cfg.trials                       = 'all' or a selection given as a 1xN vector (default = 'all')
%  cfg.(cfg.method).noiseEstimate   = 'non-parametric' or 'parametric', estimate noise
%                                       using parametric or non-parametric (default) method
%  cfg.(cfg.method).tapsmofrq       = value, smoothing parameter of noise
%                                       estimate (default = 10)
%
%
%  METHOD SPECIFIC OPTIONS AND DESCRIPTIONS
%
% ASEO
%  ASEO iteratively models single-trial event-related activity and
%  ongoing activity and gives an estimate of single-trial latency shift and
%  amplitude scaling of event-related components.
%  cfg.unit                 = 'ms' or 'sample', specify jitter and comlatency
%                              in ms or in samples (default = 'ms')
%  cfg.aseo.jitter          = value, time jitter in initial timewindow estimate
%  cfg.aseo.numiteration    = value, number of iteration (default = 1)
%  cfg.aseo.searchWindowSet = 2xN matrix, initial set of latencies of event-
%                             related components, give as [comp1start, comp1end;
%                             comp2start, comp2end] (default not specified)
%  OR
%  cfg.aseo.initcomp        = vector, initial estimate of the waveform component
%  % NOTE: multiple channels is currently not supported. if jitter is
%  specified once, jitter{k} where k is the channelnumber, is empty and
%  fails.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also <<give a list of function names, all in capitals>>

% Here come the Copyrights
%
% Here comes the Revision tag, which is auto-updated by the version control system
% $Id$


% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults                   % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init              % this will reset ft_warning and show the function help if nargin==0 and return an error
ft_preamble debug             % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar    data   % this reads the input data in case the user specified the cfg.inputfile option
ft_preamble provenance data   % this records the time and memory usage at the beginning of the function
ft_preamble trackconfig       % this converts the cfg structure in a config object, which tracks the cfg options that are being used

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% ensure that the required options are present
cfg         = ft_checkconfig(cfg, 'required', {'method'});
cfg.trials  = ft_getopt(cfg, 'trials',  'all', 1); % all trials as default
cfg.channel = ft_getopt(cfg, 'channel', 'all');

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'method', 'char', {'aseo'});

% get the method specific options
switch cfg.method
  case 'aseo'
    if ~isfield(cfg, 'aseo'), cfg.aseo = []; end
    cfg.(cfg.method).thresholdAmpH = ft_getopt(cfg.(cfg.method), 'thresholdAmpH', 0.5);
    cfg.(cfg.method).thresholdAmpL = ft_getopt(cfg.(cfg.method), 'thresholdAmpL', 0.1);
    cfg.(cfg.method).thresholdCorr = ft_getopt(cfg.(cfg.method), 'thresholdCorr', 0.05);
    cfg.(cfg.method).maxOrderAR    = ft_getopt(cfg.(cfg.method), 'maxOrderAR', 5);
    cfg.(cfg.method).noiseEstimate = ft_getopt(cfg.(cfg.method), 'noiseEstimate', 'non-parametric');
    cfg.(cfg.method).tapsmofrq     = ft_getopt(cfg.(cfg.method), 'tapsmofrq', 10);
  otherwise
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select trials of interest
tmpcfg         = [];
tmpcfg.trials  = cfg.trials;
tmpcfg.channel = cfg.channel;
data           = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% some proper error handling
if isfield(data, 'trial') && numel(data.trial)==0
  error('no trials were selected'); % this does not apply for MVAR data
end

if numel(data.label)==0
  error('no channels were selected');
end

% define general variables -> this needs to be cleaned up
fsample = data.fsample; % Sampling Frequency in Hz
period  = 1000/fsample; % Sampling Period in millsecond
nchan   = numel(data.label);
nsmp    = numel(data.time{1}); %FIXME ASSUMING FIXED TIME AXIS ACROSS ALL TRIALS

cfg.(cfg.method).fsample = fsample;
cfg.(cfg.method).searchGrid = period;
cfg.(cfg.method).sampPeri = period;
cfg.(cfg.method).nchan = nchan;
cfg.(cfg.method).nsmp = nsmp;
cfg.(cfg.method).ntrl = size(data.trial,2);


% switch over method and do some of the method specfic checks and defaulting
switch cfg.method
  
  case 'aseo'
    
    % set defaults
    waveformInitSet  = ft_getopt(cfg.aseo, 'waveformInitSet', {});
    unit             = ft_getopt(cfg.aseo, 'unit', 'ms');
    jitter           = ft_getopt(cfg.aseo, 'jitter', 50);
    numiteration     = ft_getopt(cfg.aseo, 'numiteration', 1);
    initcomp         = ft_getopt(cfg.aseo, 'initcomp', {});
    
    reconstructed = data;
    residual      = data;
    if isfield(data, 'cfg')
      reconstructed = rmfield(reconstructed, 'cfg');
      residual      = rmfield(residual, 'cfg');
    end
    for k = 1:numel(data.trial)
      reconstructed.trial{k}(:) = nan;
      residual.trial{k}(:)      = nan;
    end
    
    if isempty(waveformInitSet) && isempty(initcomp)
      error('you should supply either an initial estimate of the waveform component, or a set of latencies');
    end
    
    % if jitter universal, put it in the right format
    if length(jitter)==1
      if strcmp(unit, 'ms')
        jitter = jitter/1000*fsample; % convert jitter from ms to sample
      elseif strcmp(unit, 's')
        jitter = jitter*fsample; % convert jitter from sec to sample
      end
      jitter = repmat([-jitter jitter], [size(waveformInitSet,1) 1]);
    elseif size(jitter)<size(waveformInitSet) | size(jitter)>size(waveformInitSet)
      error('please specify cfg.aseo.jitter as a universal single value or as a matrix with size(cfg.aseo.searchWindowSet)')
    end
    cfg.aseo.jitter = jitter;
    
    if ~isempty(waveformInitSet)
      waveformInitSet = waveformInitSet';
      waveformInitSet = waveformInitSet(:);
      if strcmp(unit, 'ms')
        for k = 1:length(waveformInitSet)
          waveformInitSet(k,1) = nearest(data.time{1}, waveformInitSet(k,1)/1000); % convert unit from msec to sample
        end
      elseif strcmp(unit, 's')
        for k = 1:length(waveformInitSet)
          waveformInitSet(k,1) = nearest(data.time{1}, waveformInitSet(k,1)); % convert unit from msec to sample
        end
      end
    end
    cfg.aseo.waveformInitSet = waveformInitSet;
    cfg.aseo.unit = 'sample';
    
    % preliminaries for the inital shape of the waveform's components
    if isempty(initcomp)
      make_init = true;
      Ncomp     = size(waveformInitSet,1)./2;
    else
      make_init = false;
      % ensure it to be a cell-array
      if ~iscell(initcomp)
        initcomp = {initcomp};
      end
      Ncomp = zeros(numel(initcomp));
      for k = 1:numel(initcomp)
        Ncomp(k,1) = size(initcomp{k},2);
      end
    end
    
    % preliminaries for the jitter
    if isempty(jitter)
      jitter = cell(size(initcomp)); % MVE: what if initcomp is not yet specified?
      for k = 1:numel(jitter)
        jitter{k} = ones(Ncomp(k),1)*[-50 50];
      end
    elseif ~iscell(jitter)
      jitter = {jitter};
    end
    cfg.aseo.jitter = jitter{1};
    
    if nchan>1 && ~make_init && numel(initcomp)==1
      error('if supplying more than one channel, the initial component waveforms should be entered as a cell-array');
    end
    
    params = struct([]);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main functionality, loop over channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch cfg.method
  
  case 'aseo'
    for k = 1:nchan
      % preprocessing data
      tmp     = cellrowselect(data.trial,k);
      chandat = cat(1,tmp{:});
      
      % baseline correction
      chandat = ft_preproc_baselinecorrect(chandat, nearest(data.time{1}, -inf), nearest(data.time{1}, 0));
      
      if make_init % create initial estimate waveform from set of latencies
        avgdat   = mean(chandat, 1);
        
        % Set the initial ERP waveforms according to the preset parameters
        ncomp  = Ncomp(k);
        initcomp{k} = zeros(nsmp, ncomp);
        for m = 1:ncomp
          begsmp = waveformInitSet(2*m-1, k);
          endsmp = waveformInitSet(2*m,   k);
          if (endsmp <= begsmp), disp('Invalid input! '); end
          if begsmp<1,           begsmp = 1;              end
          if endsmp>nsmp,        endsmp = nsmp;           end
          
          tmp = avgdat(begsmp:endsmp)';
          initcomp{k}(begsmp:endsmp, m) = tmp;
        end
        initcomp{k} = initcomp{k} - repmat(mean(initcomp{k}),nsmp,1);
      end
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% ASEO algorithm
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Do zero-padding and FFT to the signal and initial waveforms
      %npad             = ceil(max(waveformInitSet(:)) -  min(waveformInitSet(:)));  % zero-padding number
      npad     = nsmp./2; % zero-padding number, I am not sure what the effect of this parameter is
      nfft     = 2.^(ceil(log2(nsmp+npad)))*2;
      initcomp_fft = fft(initcomp{k}, nfft); % Fourier transform of the initial waveform
      chandat_fft  = fft(chandat', nfft);    % Fourier transform of the signal
      output   = ft_singletrialanalysis_aseo(cfg, chandat_fft, initcomp_fft);
      
      params(k).latency    = output.lat_est;
      params(k).amplitude  = output.amp_est;
      params(k).components = output.erp_est;
      params(k).rejectflag = output.rejectflag;
      
      for m = 1:numel(data.trial)
        if output.rejectflag(m)==0
          reconstructed.trial{m}(k,:) = data.trial{m}(k,:)-output.residual(:,m)';
          residual.trial{m}(k,:)      = output.residual(:,m)';
        end
      end
    end
    reconstructed.params = params;
end


% this might involve more active checking of whether the input options
% are consistent with the data and with each other

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function

% the ft_postamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_postamble debug               % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig         % this converts the config object back into a struct and can report on the unused fields
ft_postamble previous   datain   % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble provenance dataout  % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and MATLAB version etc. to the output cfg
ft_postamble history    dataout  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar    dataout  % this saves the output data structure to disk in case the user specified the cfg.outputfile option

