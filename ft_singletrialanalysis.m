function [dataout] = ft_singletrialanalysis(cfg, data)

% FT_SINGLETRIALANALYSIS computes a model of single trial event-
% related activity
%
% Use as
%   [dataout] = ft_singletrialanalysis(cfg, data)
% where data is single-channel raw data as obtained by FT_PREPROCESSING
% and cfg is a configuration structure according to
%
%  cfg.method  = method for estimating event-related activity
%                 'aseo', analysis of single-trial ERP and ongoing
%                         activity (according to Xu et al, 2009)
%                 'gbve', graph-based variability estimation
%                         (according to Gramfort et al, IEEE TBME 2009)
%  cfg.channel = Nx1 cell-array with selection of channels (default = 'all'),
%                see FT_CHANNELSELECTION for details
%  cfg.trials  = 'all' or a selection given as a 1xN vector (default = 'all')
%  cfg.output  = 'model', or 'residual', which returns the modelled data,
%                or the residuals.
%
% METHOD SPECIFIC OPTIONS AND DESCRIPTIONS
%
% ASEO
%  ASEO iteratively models single-trial event-related activity and
%  ongoing activity and gives an estimate of single-trial latency shift and
%  amplitude scaling of event-related components. It has the following options:
%
%  cfg.aseo.noiseEstimate   = 'non-parametric' or 'parametric', estimate noise
%                               using parametric or non-parametric (default) method
%  cfg.aseo.tapsmofrq       = value, smoothing parameter of noise
%                               estimate (default = 10)
%  cfg.aseo.jitter          = value, time jitter in initial timewindow
%                              estimate (in seconds). default 0.050 seconds
%  cfg.aseo.numiteration    = value, number of iteration (default = 1)
%  cfg.aseo.waveformInitSet = 2xN matrix, initial set of latencies in seconds of event-
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
% See also FT_SINGLETRIALANALYSIS_ASEO

% Copyright (C) 2018, Jan-Mathijs Schoffelen DCCN
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure that the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% ensure that the required options are present
cfg         = ft_checkconfig(cfg, 'required', {'method'});
cfg.trials  = ft_getopt(cfg, 'trials',  'all', 1); % all trials as default
cfg.channel = ft_getopt(cfg, 'channel', 'all');
cfg.output  = ft_getopt(cfg, 'output',  'model');
% ensure that the options are valid
cfg = ft_checkopt(cfg, 'method', 'char', {'aseo' 'gbve'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select trials of interest
tmpcfg = keepfields(cfg, {'trials' 'channel' 'showcallinfo'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% some error checks
if isfield(data, 'trial') && numel(data.trial)==0, ft_error('no trials were selected'); end
if numel(data.label)==0, ft_error('no channels were selected'); end

switch cfg.method  
  case 'aseo'
    % define general variables that are used locally
    fsample = data.fsample; % Sampling Frequency in Hz
    nchan   = numel(data.label);
    nsmp    = numel(data.time{1}); %FIXME ASSUMING FIXED TIME AXIS ACROSS ALL TRIALS

    % setting a bunch of options, to be passed on to the lower level function
    if ~isfield(cfg, 'aseo'), cfg.aseo = []; end 
    cfg.aseo.thresholdAmpH = ft_getopt(cfg.aseo, 'thresholdAmpH', 0.5);
    cfg.aseo.thresholdAmpL = ft_getopt(cfg.aseo, 'thresholdAmpL', 0.1);
    cfg.aseo.thresholdCorr = ft_getopt(cfg.aseo, 'thresholdCorr', 0.05);
    cfg.aseo.maxOrderAR    = ft_getopt(cfg.aseo, 'maxOrderAR',    5);
    cfg.aseo.noiseEstimate = ft_getopt(cfg.aseo, 'noiseEstimate', 'non-parametric');
    cfg.aseo.numiteration  = ft_getopt(cfg.aseo, 'numiteration',  1);
    cfg.aseo.tapsmofrq     = ft_getopt(cfg.aseo, 'tapsmofrq',     10);
    cfg.aseo.fsample       = fsample;
    cfg.aseo.searchGrid    = 1000./fsample; % defined in ms
    cfg.aseo.sampPeri      = 1000./fsapmle; % defined in ms;
    cfg.aseo.nchan         = nchan;
    cfg.aseo.nsmp          = nsmp;
    cfg.aseo.ntrl          = numel(data.trial);
        
    % deal with the different ways with which the initial waveforms can be defined
    waveformInitSet  = ft_getopt(cfg.aseo, 'waveformInitSet', {});
    jitter           = ft_getopt(cfg.aseo, 'jitter', 0.020);
    initcomp         = ft_getopt(cfg.aseo, 'initcomp', {});
        
    if isempty(waveformInitSet) && isempty(initcomp)
      ft_error('you should supply either an initial estimate of the waveform component, or a set of latencies');
    end
 
    % if jitter universal, put it in the right format
    if length(jitter)==1
      jitter = jitter*fsample; % convert jitter from sec to sample -> this is inconsistent with the definition of another time parameter in ms
      jitter = repmat([-jitter jitter], [size(waveformInitSet,1) 1]);
    elseif size(jitter)==size(waveformInitSet)
      jitter = jitter*fsample;
    elseif size(jitter)<size(waveformInitSet) || size(jitter)>size(waveformInitSet)
      ft_error('please specify cfg.aseo.jitter as a universal single value or as a matrix with size(cfg.aseo.searchWindowSet)')
    end
    cfg.aseo.jitter = jitter;
        
    if ~isempty(waveformInitSet)
      waveformInitSet = waveformInitSet';
      waveformInitSet = waveformInitSet(:);
      % convert seconds into samples
      for k = 1:length(waveformInitSet)
        waveformInitSet(k,1) = nearest(data.time{1}, waveformInitSet(k,1)); % convert unit from sec to sample
      end
    end
    cfg.aseo.waveformInitSet = waveformInitSet;
    cfg.aseo.unit            = 'sample';
        
    % preliminaries for the inital shape of the waveform's components
    if isempty(initcomp)
      % this results in the initial components to be estimated from the ERP, given the specification of 
      % approximate latencies
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
        jitter{k} = ones(Ncomp(k),1)*[-cfg.aseo.jitter cfg.aseo.jitter];
      end
    elseif ~iscell(jitter)
      jitter = {jitter};
    end
    cfg.aseo.jitter = jitter{1};
        
    if nchan>1 && ~make_init && numel(initcomp)==1
      ft_error('if supplying more than one channel, the initial component waveforms should be entered as a cell-array');
    end
        
    % initialize the output data
    dataout = removefields(data, 'cfg');
    for k = 1:numel(data.trial)
      dataout.trial{k}(:) = nan;
    end
        
    % initialize the struct that will contain the output parameters
    params = struct([]);
        
    for k = 1:nchan
      % preprocessing data
      tmp     = cellrowselect(data.trial,k);
      chandat = cat(1,tmp{:});
            
      % baseline correction
      chandat = ft_preproc_baselinecorrect(chandat, nearest(data.time{1}, -inf), nearest(data.time{1}, 0));
      
      % create initial estimate waveform from set of latencies
      if make_init
        avgdat   = mean(chandat, 1);
                
        % Set the initial ERP waveforms according to the preset parameters
        ncomp       = Ncomp(k);
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
            
      % do zero-padding and FFT to the signal and initial waveforms
      %npad             = ceil(max(waveformInitSet(:)) -  min(waveformInitSet(:)));  % zero-padding number
      npad     = nsmp./2; % zero-padding number, I am not sure what the effect of this parameter is
      nfft     = 2.^(ceil(log2(nsmp+npad)))*2;
      initcomp_fft = fft(initcomp{k}, nfft); % Fourier transform of the initial waveform
      chandat_fft  = fft(chandat', nfft);    % Fourier transform of the signal
      output   = ft_singletrialanalysis_aseo(cfg, chandat_fft, initcomp_fft);
      
      params(k).latency    = output.lat_est/fsample;
      params(k).amplitude  = output.amp_est_unscaled;
      params(k).components = output.erp_est;
      params(k).rejectflag = output.rejectflag;
      
      for m = 1:numel(data.trial)
        if output.rejectflag(m)==0
          switch cfg.output
            case 'model'
              dataout.trial{m}(k,:) = data.trial{m}(k,:)-output.residual(:,m)';
            case 'residual'
              dataout.trial{m}(k,:) = output.residual(:,m)';
          end
        end
      end
    end
    
case 'gbve'
  ft_hastoolbox('lagextraction', 1);
  ft_hastoolbox('eeglab',        1); % because the low-level code might use a specific moving average function from EEGLAB
  ft_hastoolbox('cellfunction',  1);
  
  if ~isfield(cfg, 'gbve'), cfg.gbve = []; end
  cfg.gbve.NORMALIZE_DATA    = ft_getopt(cfg.gbve, 'NORMALIZE_DATA',     true);
  cfg.gbve.CENTER_DATA       = ft_getopt(cfg.gbve, 'CENTER_DATA',        false);
  cfg.gbve.USE_ADAPTIVE_SIGMA= ft_getopt(cfg.gbve, 'USE_ADAPTIVE_SIGMA', false);
  cfg.gbve.sigma             = ft_getopt(cfg.gbve, 'sigma',    0.01:0.01:0.2);
  cfg.gbve.distance          = ft_getopt(cfg.gbve, 'distance', 'corr2');
  cfg.gbve.alpha             = ft_getopt(cfg.gbve, 'alpha',    [0 0.001 0.01 0.1]);
  cfg.gbve.exponent          = ft_getopt(cfg.gbve, 'exponent', 1);
  cfg.gbve.use_maximum       = ft_getopt(cfg.gbve, 'use_maximum', 1); % consider the positive going peak
  cfg.gbve.show_pca          = ft_getopt(cfg.gbve, 'show_pca',          false);
  cfg.gbve.show_trial_number = ft_getopt(cfg.gbve, 'show_trial_number', false);
  cfg.gbve.verbose           = ft_getopt(cfg.gbve, 'verbose',           true);
  cfg.gbve.disp_log          = ft_getopt(cfg.gbve, 'disp_log',          false);
  cfg.gbve.latency           = ft_getopt(cfg.gbve, 'latency',  [-inf inf]);
  cfg.gbve.xwin              = ft_getopt(cfg.gbve, 'xwin',     1); % default is a bit of smoothing
  
  nchan = numel(data.label);
  ntrl  = numel(data.trial);
  
  tmin  = nearest(data.time{1}, cfg.gbve.latency(1));
  tmax  = nearest(data.time{1}, cfg.gbve.latency(2));

  % initialize the struct that will contain the output parameters
  dataout = removefields(data, 'cfg');
  params  = struct([]);
  options = cfg.gbve;
  for k = 1:nchan
    % preprocessing data
    tmp     = cellrowselect(data.trial,k);
    chandat = cat(1,tmp{:});
    points  = chandat(:,tmin:tmax);
    
    % perform a loop across alpha values, cross validation
    alphas = options.alpha;

    if length(alphas) > 1 % Use Cross validation error if multiple alphas are specified
      best_CVerr = -Inf;

      K = 5;
      disp(['--- Running K Cross Validation (K = ',num2str(K),')']);

      block_idx = fix(linspace(1, ntrl, K+1)); % K cross validation
      for jj=1:length(alphas)
        options.alpha = alphas(jj);

        CVerr = 0;
        for kk = 1:K
          bidx = block_idx(jj):block_idx(jj+1);
          idx = 1:ntrl;
          idx(bidx) = [];

          data_k       = chandat(idx,:);
          points_k     = points(idx,:);
          [order,lags] = extractlag(points_k,options);

          data_reordered = data_k(order,:);
          lags           = lags + tmin;
          [data_aligned, ~] = perform_realign(data_reordered, data.time{1}, lags);
          data_aligned(~isfinite(data_aligned)) = nan;
          ep_evoked = nanmean(data_aligned);
          ep_evoked = ep_evoked ./ norm(ep_evoked);

          data_k = chandat(bidx,:);
          data_norm = sqrt(sum(data_k.^2,2));
          data_k = diag(1./data_norm)*data_k;
          data_k(data_norm==0,:) = 0;
          
          for pp=1:length(bidx)
            c     = xcorr(ep_evoked,data_k(pp,:));
            CVerr = CVerr + max(c(:));
          end
        end

        CVerr = CVerr/ntrl;

        if CVerr > best_CVerr
          best_CVerr = CVerr;
          best_alpha = alphas(jj);
        end
      end
      options.alpha = best_alpha;
    end

    if options.use_maximum
      [order,lags] = extractlag( points, options );
    else
      [order,lags] = extractlag( -points, options );
    end
    disp(['---------- Using alpha = ',num2str(options.alpha)]);
    data_reordered = chandat(order,:);
    lags = lags + tmin;
    [data_aligned] = perform_realign(data_reordered, data.time{1}, lags );
    data_aligned(~isfinite(data_aligned)) = nan;
    
    [~,order_inv] = sort(order);
    lags_no_order = lags(order_inv);
    data_aligned  = data_aligned(order_inv,:);
      
    params(k).lags = [lags_no_order data.time{1}(lags_no_order)'];
    switch cfg.output
      case 'model'
        tmp = mat2cell(data_aligned, ones(1,size(data_aligned,1)), size(data_aligned,2))';
        dataout.trial = cellrowassign(dataout.trial, tmp, k);
      case 'residual'
        % to be done
        error('not yet implemented');
    end
    
    

  end
end

dataout.params = params;
dataout.cfg    = cfg;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

