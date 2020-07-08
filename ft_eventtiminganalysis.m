function [dataout] = ft_eventtiminganalysis(cfg, data)

% FT_EVENTTIMINGANALYSIS computes a model of single trial event-
% related activity, by estimating per trial the latency (and
% amplitude) of event-related signal components.
%
% Use as
%   [dataout] = ft_eventtiminganalysis(cfg, data)
% where data is single-channel raw data as obtained by FT_PREPROCESSING
% and cfg is a configuration structure according to
%
%   cfg.method  = method for estimating event-related activity
%                 'aseo', analysis of single-trial ERP and ongoing
%                         activity (according to Xu et al, 2009)
%                 'gbve', graph-based variability estimation
%                         (according to Gramfort et al, IEEE TBME 2009)
%   cfg.channel = Nx1 cell-array with selection of channels (default = 'all'),
%                 see FT_CHANNELSELECTION for details
%   cfg.trials  = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.output  = 'model', or 'residual', which returns the modelled data,
%                 or the residuals.
%
% Method specific options are specified in the appropriate substructure.
%
% For the ASEO method, the following options can be specified:
%   cfg.aseo.noiseEstimate   = 'non-parametric' or 'parametric', estimate noise
%                              using parametric or non-parametric (default) method
%   cfg.aseo.tapsmofrq       = value, smoothing parameter of noise for
%                              nonparametric estimation (default = 5)
%   cfg.aseo.jitter          = value, time jitter in initial timewindow
%                              estimate (in seconds). default 0.050 seconds
%   cfg.aseo.numiteration    = value, number of iteration (default = 1)
%   cfg.aseo.initlatency     = Nx2 matrix, initial set of latencies in seconds of event-
%                              related components, give as [comp1start, comp1end;
%                              comp2start, comp2end] (default not
%                              specified). For multiple channels it should
%                              be a cell-array, one matrix per channel
%  Alternatively, rather than specifying a (set of latencies), one can also
%  specify:
%
%   cfg.aseo.initcomp        = vector, initial estimate of the waveform
%                              components. For multiple channels it should
%                              be a cell-array, one matrix per channel.
%
% For the GBVE method, the following options can be specified:
%   cfg.gbve.sigma             = vector, range of sigma values to explore in
%                                cross-validation loop (default: 0.01:0.01:0.2)
%   cfg.gbve.distance          = scalar, distance metric to use as
%                                evaluation criterion, see plugin code for
%                                more informatoin
%   cfg.gbve.alpha             = vector, range of alpha values to explor in
%                                cross-validation loop (default: [0 0.001 0.01 0.1])
%   cfg.gbve.exponent          = scalar, see plugin code for information
%   cfg.gbve.use_maximum       = boolean, (default: 1) consider the positive going peak
%   cfg.gbve.show_pca          = boolean, see plugin code (default 0)
%   cfg.gbve.show_trial_number = boolean, see plugin code (default 0)
%   cfg.gbve.verbose           = boolean (default: 1)
%   cfg.gbve.disp_log          = boolean, see plugin code (default 0)
%   cfg.gbve.latency           = vector [min max], latency range in s
%                                (default: [-inf inf])
%   cfg.gbve.xwin              = scalar smoothing parameter for moving
%                                average smoothing (default: 1), see
%                                eeglab's movav function for more
%                                information.
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

% Copyright (C) 2018-2019, Jan-Mathijs Schoffelen DCCN
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
    nsample = numel(data.time{1}); %FIXME ASSUMING FIXED TIME AXIS ACROSS ALL TRIALS

    % setting a bunch of options, to be passed on to the lower level function
    if ~isfield(cfg, 'aseo'), cfg.aseo = []; end
    cfg.aseo.thresholdAmpH = ft_getopt(cfg.aseo, 'thresholdAmpH', 0.5);
    cfg.aseo.thresholdAmpL = ft_getopt(cfg.aseo, 'thresholdAmpL', 0.1);
    cfg.aseo.thresholdCorr = ft_getopt(cfg.aseo, 'thresholdCorr', 0.2);
    cfg.aseo.maxOrderAR    = ft_getopt(cfg.aseo, 'maxOrderAR',    5);
    cfg.aseo.noiseEstimate = ft_getopt(cfg.aseo, 'noiseEstimate', 'nonparametric');
    cfg.aseo.numiteration  = ft_getopt(cfg.aseo, 'numiteration',  1);
    cfg.aseo.tapsmofrq     = ft_getopt(cfg.aseo, 'tapsmofrq',     5);
    cfg.aseo.fsample       = fsample;
    cfg.aseo.nsample       = nsample;
    cfg.aseo.pad           = ft_getopt(cfg.aseo, 'pad', (2.*nsample)/fsample);
    
    % deal with the different ways with which the initial waveforms can be defined
    initlatency      = ft_getopt(cfg.aseo, 'initlatency', {});
    initcomp         = ft_getopt(cfg.aseo, 'initcomp',    {});
    jitter           = ft_getopt(cfg.aseo, 'jitter',      0.050); % half temporal width of shift in s
        
    if isempty(initlatency) && isempty(initcomp)
      ft_error('for the ASEO method you should supply either an initial estimate of the waveform component, or a set of latencies');
    elseif ~isempty(initlatency)
      % this takes precedence, and should contain per channel the begin and
      % end points of the subwindows in time, based on which the initial
      % subcomponents are estimated
    
      % ensure it to be a cell-array if the input is a matrix
      if ~iscell(initlatency)
        initlatency = repmat({initlatency},[1 nchan]);
      end
      make_init = true;
    elseif ~isempty(initcomp)
      % ensure it to be a cell-array if the input is a matrix
      if ~iscell(initcomp)
        initcomp = repmat({initcomp}, [1 nchan]);
      end
      make_init = false;
    end
    
    if make_init
      assert(numel(initlatency)==nchan);
      for k = 1:nchan
        % preprocessing data
        tmp     = cellrowselect(data.trial,k);
        chandat = cat(1,tmp{:});
        chandat = ft_preproc_baselinecorrect(chandat, nearest(data.time{1}, -inf), nearest(data.time{1}, 0));
        avgdat  = nanmean(chandat, 1);
                
        % set the initial ERP waveforms according to the preset parameters
        ncomp       = size(initlatency{k},1);
        initcomp{k} = zeros(nsample, ncomp);
        for m = 1:ncomp
          begsmp = nearest(data.time{1},initlatency{k}(m, 1));
          endsmp = nearest(data.time{1},initlatency{k}(m, 2));
          if begsmp<1,       begsmp = 1;       end
          if endsmp>nsample, endsmp = nsample; end
                 
          tmp = avgdat(begsmp:endsmp)';
          initcomp{k}(begsmp:endsmp, m) = tmp;
        end
        initcomp{k} = initcomp{k} - repmat(mean(initcomp{k}),nsample,1);
      end
    else
      assert(numel(initcomp)==nchan);
    end
    
    if ~iscell(jitter)
      jitter = repmat({jitter}, [1 nchan]);
    end
    
    for k = 1:numel(jitter)
      if ~isempty(jitter{k})
        if size(jitter{k},1)~=size(initcomp{k},2), jitter{k} = repmat(jitter{k}(1,:),[size(initcomp{k},2) 1]); end
      end
    end
    
    % initialize the output data
    dataout = removefields(data, 'cfg');
    for k = 1:numel(data.trial)
      dataout.trial{k}(:) = nan;
    end
        
    % initialize the struct that will contain the output parameters
    params = struct([]);
    
    % do the actual computations
    for k = 1:nchan
      % preprocessing data
      tmp     = cellrowselect(data.trial,k);
      chandat = cat(1,tmp{:});
            
      % baseline correction
      chandat = ft_preproc_baselinecorrect(chandat, nearest(data.time{1}, -inf), nearest(data.time{1}, 0));
      
      % do zero-padding and FFT to the signal and initial waveforms
      npad         = cfg.aseo.pad*fsample;    % length of data + zero-padding number
      nfft         = 2.^(ceil(log2(npad)))*2;
      initcomp_fft = fft(initcomp{k}, nfft);  % Fourier transform of the initial waveform
      chandat_fft  = fft(chandat', nfft);     % Fourier transform of the signal
      
      cfg.aseo.jitter = jitter{k};
      output       = ft_singletrialanalysis_aseo(cfg, chandat_fft, initcomp_fft);
      
      params(k).latency    = output(end).lat_est./fsample;
      params(k).amplitude  = output(end).amp_est;
      params(k).components = output(end).erp_est;
      params(k).rejectflag = output(end).rejectflag;
      params(k).noise      = output(end).noise;
      
      for m = 1:numel(data.trial)
        if output(end).rejectflag(m)==0
          switch cfg.output
            case 'model'
              dataout.trial{m}(k,:) = data.trial{m}(k,:)-output(end).residual(:,m)';
            case 'residual'
              dataout.trial{m}(k,:) = output(end).residual(:,m)';
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
  cfg.gbve.nfold             = ft_getopt(cfg.gbve, 'nfold',    5);
  
  nchan = numel(data.label);
  ntrl  = numel(data.trial);
  
  tmin  = nearest(data.time{1}, cfg.gbve.latency(1));
  tmax  = nearest(data.time{1}, cfg.gbve.latency(2));

  % initialize the struct that will contain the output parameters
  dataout = removefields(data, 'cfg');
  params  = struct([]);
  for k = 1:nchan
    % preprocessing data
    options = cfg.gbve;
  
    fprintf('--- Processing channel %d\n',k);
    
    tmp     = cellrowselect(data.trial,k);
    chandat = cat(1,tmp{:});
    points  = chandat(:,tmin:tmax);
    
    % perform a loop across alpha values, cross validation
    alphas = options.alpha;
    
    if length(alphas) > 1 % Use Cross validation error if multiple alphas are specified
      best_CVerr = -Inf;

      K = cfg.gbve.nfold;
      disp(['--- Running K Cross Validation (K = ',num2str(K),')']);

      block_idx = fix(linspace(1, ntrl, K+1)); % K cross validation
      for jj=1:length(alphas)
        options.alpha = alphas(jj);

        CVerr = 0;
        for kk = 1:K
          bidx = block_idx(kk):block_idx(kk+1);
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
      
    params(k).latency = data.time{1}(lags_no_order)';
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
