function [mvardata] = ft_mvaranalysis(cfg, data)

% FT_MVARANALYSIS performs multivariate autoregressive modeling on
% time series data over multiple trials.
%
% Use as
%   [mvardata] = ft_mvaranalysis(cfg, data)
%
% The input data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration depends on the type
% of computation that you want to perform.
% The output is a data structure of datatype 'mvar' which contains the
% multivariate autoregressive coefficients in the field coeffs, and the
% covariance of the residuals in the field noisecov.
%
% The configuration should contain:
%   cfg.method     = the name of the toolbox containing the function for the
%                     actual computation of the ar-coefficients
%                     this can be 'biosig' (default) or 'bsmart'
%                     you should have a copy of the specified toolbox in order
%                     to use mvaranalysis (both can be downloaded directly).
%   cfg.mvarmethod = scalar (only required when cfg.method = 'biosig').
%                     default is 2, relates to the algorithm used for the
%                     computation of the AR-coefficients by mvar.m
%   cfg.order      = scalar, order of the autoregressive model (default=10)
%   cfg.channel    = 'all' (default) or list of channels for which an mvar model
%                     is fitted. (Do NOT specify if cfg.channelcmb is
%                     defined)
%   cfg.channelcmb = specify channel combinations as a
%                     two-column cell-array with channels in each column between
%                     which a bivariate model will be fit (overrides
%                     cfg.channel)
%   cfg.keeptrials = 'no' (default) or 'yes' specifies whether the coefficients
%                     are estimated for each trial separately, or on the
%                     concatenated data
%   cfg.jackknife  = 'no' (default) or 'yes' specifies whether the coefficients
%                     are estimated for all leave-one-out sets of trials
%   cfg.zscore     = 'no' (default) or 'yes' specifies whether the channel data
%                      are z-transformed prior to the model fit. This may be
%                      necessary if the magnitude of the signals is very different
%                      e.g. when fitting a model to combined MEG/EMG data
%   cfg.demean     = 'yes' (default) or 'no' explicit removal of DC-offset
%   cfg.ems        = 'no' (default) or 'yes' explicit removal ensemble mean
%
% ft_mvaranalysis can be used to obtain one set of coefficients across
% all time points in the data, also when the trials are of varying length.
%
% ft_mvaranalysis can be also used to obtain time-dependent sets of
% coefficients based on a sliding window. In this case the input cfg
% should contain:
%
%   cfg.t_ftimwin = the width of the sliding window on which the coefficients
%                    are estimated
%   cfg.toi       = [t1 t2 ... tx] the time points at which the windows are
%                    centered
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_PREPROCESSING, FT_SOURCESTATISTICS, FT_FREQSTATISTICS,
% FT_TIMELOCKSTATISTICS

% Undocumented local options:
%   cfg.keeptapers
%   cfg.taper
%   cfg.output      = 'parameters', 'model', 'residual'. If 'parameters' is
%     specified, the output is a mdata data structure, containing the
%     coefficients and the noise covariance. If 'model' or 'residual' is
%     specified, the output is a data structure containing either the
%     modeled time series, or the residuals. This is only supported when
%     the model is estimated across the whole time range.

% Copyright (C) 2009, Jan-Mathijs Schoffelen
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% this must be done prior to "ft_preamble init" which merges the cfg with the global ft_default
if isfield(cfg, 'toolbox') && any(strcmp(cfg.toolbox, {'bsmart', 'biosig'}))
  ft_warning('please use cfg.method instead of cfg.toolbox');
  % cfg.toolbox is used in ft_default
  cfg.method = cfg.toolbox;
  cfg = rmfield(cfg, 'toolbox');
end

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'hassampleinfo', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed', {'blcwindow', 'baselinewindow'});

% set default configuration options
cfg.method     = ft_getopt(cfg, 'method', 'biosig');
cfg.mvarmethod = ft_getopt(cfg, 'mvarmethod', 2); % only relevant for biosig
cfg.order      = ft_getopt(cfg, 'order',      10);
cfg.channel    = ft_getopt(cfg, 'channel',    'all');
cfg.keeptrials = ft_getopt(cfg, 'keeptrials', 'no');
cfg.jackknife  = ft_getopt(cfg, 'jackknife',  'no');
cfg.zscore     = ft_getopt(cfg, 'zscore',     'no');
cfg.feedback   = ft_getopt(cfg, 'feedback',   'textbar');
cfg.demean     = ft_getopt(cfg, 'demean',     'yes');
cfg.ems        = ft_getopt(cfg, 'ems',        'no');
cfg.toi        = ft_getopt(cfg, 'toi',        []);
cfg.t_ftimwin  = ft_getopt(cfg, 't_ftimwin',  []);
cfg.keeptapers = ft_getopt(cfg, 'keeptapers', 'yes');
cfg.taper      = ft_getopt(cfg, 'taper',      'rectwin');
cfg.univariate = ft_getopt(cfg, 'univariate', 0);
cfg.output     = ft_getopt(cfg, 'output',     'parameters');

% check that cfg.channel and cfg.channelcmb are not both specified
if ~any(strcmp(cfg.channel, 'all')) && isfield(cfg, 'channelcmb')
  ft_warning('cfg.channelcmb defined, overriding cfg.channel setting and computing over bivariate pairs');
else
  % select trials of interest
  tmpcfg = [];
  tmpcfg.channel = cfg.channel;
  data = ft_selectdata(tmpcfg, data);
  % restore the provenance information
  [cfg, data] = rollback_provenance(cfg, data);
end

% check whether the requested toolbox is present and check the configuration
switch cfg.method
  case 'biosig'
    % check the configuration
    cfg = ft_checkconfig(cfg, 'required', 'mvarmethod');
    ft_hastoolbox('biosig', 1);
    nnans = cfg.order;
  case 'bsmart'
    ft_hastoolbox('bsmart', 1);
    nnans = 0;
  otherwise
    ft_error('toolbox %s is not yet supported', cfg.method);
end

if isempty(cfg.toi) && isempty(cfg.t_ftimwin)
  % fit model to entire data segment
  
  % check whether this is allowed
  nsmp = cellfun('size', data.trial, 2);
  if all(nsmp==nsmp(1))
    oktoolbox = {'bsmart' 'biosig'};
  else
    oktoolbox = 'biosig'; % bsmart does not work with variable trials
  end
  
  if ~ismember(cfg.method, oktoolbox)
    error('fitting the mvar-model is not possible with the ''%s'' toolbox',cfg.method);
  end
  latency = [-inf inf];
elseif ~isempty(cfg.toi) && ~isempty(cfg.t_ftimwin)
  % do sliding window approach
  for k = 1:numel(cfg.toi)
    latency(k,:) = cfg.toi(k) + cfg.t_ftimwin.*[-0.5 0.5] + [0 -1./data.fsample];
  end
else
  ft_error('cfg should contain both cfg.toi and cfg.t_ftimwin');
end


keeprpt  = istrue(cfg.keeptrials);
keeptap  = istrue(cfg.keeptapers);
dojack   = istrue(cfg.jackknife);
dozscore = istrue(cfg.zscore);
dobvar   = isfield(cfg,           'channelcmb');
dounivariate = istrue(cfg. univariate);

if ~keeptap, ft_error('not keeping tapers is not possible yet'); end
if dojack && keeprpt, ft_error('you cannot simultaneously keep trials and do jackknifing'); end

ntrl     = length(data.trial);
ntoi     = size(latency, 1);

if ~dobvar
  chanindx = match_str(data.label, cfg.channel);
  nchan    = length(chanindx);
  label    = data.label(chanindx);
  
  ncmb     = nchan*nchan;
  cmbindx1 = repmat(chanindx(:),  [1 nchan]);
  cmbindx2 = repmat(chanindx(:)', [nchan 1]);
  labelcmb = [data.label(cmbindx1(:)) data.label(cmbindx2(:))];
else
  cfg.channelcmb = ft_channelcombination(cfg.channelcmb, data.label);
  cmbindx        = zeros(size(cfg.channelcmb));
  for k = 1:size(cmbindx,1)
    [tmp, cmbindx(k,:)] = match_str(cfg.channelcmb(k,:)', data.label);
  end
  
  
  nchan    = 2;
  label    = data.label(cmbindx);
  
  ncmb     = nchan*nchan;
  labelcmb = cell(0,2);
  cmb      = cfg.channelcmb;
  
  for k = 1:size(cmbindx,1)
    labelcmb{end+1,1} = [cmb{k,1},'[',cmb{k,1},cmb{k,2},']'];
    labelcmb{end  ,2} = [cmb{k,1},'[',cmb{k,1},cmb{k,2},']'];
    labelcmb{end+1,1} = [cmb{k,2},'[',cmb{k,1},cmb{k,2},']'];
    labelcmb{end  ,2} = [cmb{k,1},'[',cmb{k,1},cmb{k,2},']'];
    labelcmb{end+1,1} = [cmb{k,1},'[',cmb{k,1},cmb{k,2},']'];
    labelcmb{end  ,2} = [cmb{k,2},'[',cmb{k,1},cmb{k,2},']'];
    labelcmb{end+1,1} = [cmb{k,2},'[',cmb{k,1},cmb{k,2},']'];
    labelcmb{end  ,2} = [cmb{k,2},'[',cmb{k,1},cmb{k,2},']'];
  end
end

tfwin   = round(cfg.t_ftimwin.*data.fsample);
%---think whether this makes sense at all
if strcmp(cfg.taper, 'dpss')
  % create a sequence of DPSS (Slepian) tapers
  % ensure that the input arguments are double precision
  tap = double_dpss(tfwin,tfwin*(cfg.tapsmofrq./data.fsample))';
  tap = tap(1,:); %only use first 'zero-order' taper
elseif strcmp(cfg.taper, 'sine')
  tap = sine_taper(tfwin, tfwin*(cfg.tapsmofrq./data.fsample))';
  tap = tap(1,:);
else
  tap = window(cfg.taper, tfwin)';
  tap = tap./norm(tap);
end
ntap = size(tap,1);


%---preprocess data if necessary -> changed 20150224, JM does not think
%this step is necessary: it creates problems downstream if the time axes of
%the trials are different
%---cut off the uninteresting data segments
%tmpcfg        = [];
%tmpcfg.toilim = cfg.toi([1 end]) + cfg.t_ftimwin.*[-0.5 0.5];
%data          = ft_redefinetrial(tmpcfg, data);

%---demean
if strcmp(cfg.demean, 'yes')
  tmpcfg           = [];
  tmpcfg.demean    = 'yes';
  tmpcfg.baselinewindow = latency([1 end]);
  data             = ft_preprocessing(tmpcfg, data);
else
  %do nothing
end

%---ensemble mean subtraction
if strcmp(cfg.ems, 'yes')
  % to be implemented
  error('ensemble mean subtraction is not yet implemented here');
end

%---zscore
if dozscore
  zwindow = latency([1 end]);
  sumval  = 0;
  sumsqr  = 0;
  numsmp  = 0;
  trlindx = [];
  for k = 1:ntrl
    begsmp = nearest(data.time{k}, zwindow(1));
    endsmp = nearest(data.time{k}, zwindow(2));
    if endsmp>=begsmp
      sumval  = sumval + sum(data.trial{k}(:, begsmp:endsmp),    2);
      sumsqr  = sumsqr + sum(data.trial{k}(:, begsmp:endsmp).^2, 2);
      numsmp  = numsmp + endsmp - begsmp + 1;
      trlindx = [trlindx; k];
    end
  end
  datavg = sumval./numsmp;
  datstd = sqrt(sumsqr./numsmp - (sumval./numsmp).^2);
  
  data.trial = data.trial(trlindx);
  data.time  = data.time(trlindx);
  ntrl       = length(trlindx);
  for k = 1:ntrl
    rvec          = ones(1,size(data.trial{k},2));
    data.trial{k} = (data.trial{k} - datavg*rvec)./(datstd*rvec);
  end
else
  %do nothing
end

%---generate time axis
maxtim = -inf;
mintim = inf;
for k = 1:ntrl
  maxtim = max(maxtim, data.time{k}(end));
  mintim = min(mintim, data.time{k}(1));
end
timeaxis = mintim:1/data.fsample:maxtim;

%---allocate memory
if dobvar && (keeprpt || dojack)
  % not yet implemented
  ft_error('doing bivariate model fits in combination with multiple replicates is not yet possible');
elseif dobvar
  coeffs   = zeros(1, 2*nchan,  size(cmbindx,1), cfg.order, ntoi, ntap);
  noisecov = zeros(1, 2*nchan,  size(cmbindx,1),            ntoi, ntap);
elseif dounivariate && (keeprpt || dojack)
  error('doing univariate model fits in combination with multiple replicates is not yet possible');
elseif dounivariate
  coeffs   = zeros(1, nchan, cfg.order, ntoi, ntap);
  noisecov = zeros(1, nchan,            ntoi, ntap);
elseif (keeprpt || dojack)
  coeffs   = zeros(length(data.trial), nchan, nchan, cfg.order, ntoi, ntap);
  noisecov = zeros(length(data.trial), nchan, nchan,            ntoi, ntap);
else
  coeffs   = zeros(1, nchan, nchan, cfg.order, ntoi, ntap);
  noisecov = zeros(1, nchan, nchan,            ntoi, ntap);
end

%---loop over the tois
ft_progress('init', cfg.feedback, 'computing AR-model');
for j = 1:ntoi
  
  if ~isequal(latency(j,:),[-inf inf])
    ft_progress(j/ntoi, 'processing timewindow %d from %d\n', j, ntoi);
    tmpcfg = [];
    tmpcfg.toilim = latency(j,:);
    tmpdata = ft_redefinetrial(tmpcfg, data);
  else
    tmpdata = data;
  end
  
  tmpnsmp = cellfun('size', tmpdata.trial, 2);
  tfwin   = tmpnsmp(1);
  %---think whether this makes sense at all
  if strcmp(cfg.taper, 'dpss')
    % create a sequence of DPSS (Slepian) tapers
    % ensure that the input arguments are double precision
    tap = double_dpss(tfwin,tfwin*(cfg.tapsmofrq./tmpdata.fsample))';
    tap = tap(1,:); %only use first 'zero-order' taper
  elseif strcmp(cfg.taper, 'sine')
    tap = sine_taper(tfwin, tfwin*(cfg.tapsmofrq./tmpdata.fsample))';
    tap = tap(1,:);
  else
    tap = window(cfg.taper, tfwin)';
    tap = tap./norm(tap);
  end
  ntap = size(tap,1);

  
  if ntoi>1 && strcmp(cfg.method, 'bsmart')
    % ensure all segments to be of equal length
    if ~all(tmpnsmp==tmpnsmp(1))
      error('the epochs are of unequal length, possibly due to numerical time axis issues, or due to partial artifacts, use cfg.method=''biosig''');
    end
    ix         = find(tmpnsmp==mode(tmpnsmp), 1, 'first');
    cfg.toi(j) = mean(tmpdata.time{ix}([1 end]))+0.5./data.fsample; %FIXME think about this
  end
  
  
  %---create cell-array indexing which original trials should go into each replicate
  rpt  = {};
  nrpt = numel(tmpdata.trial);
  if dojack
    rpt = cell(nrpt,1);
    for k = 1:nrpt
      rpt{k,1} = setdiff(1:nrpt,k);
    end
  elseif keeprpt
    for k = 1:nrpt
      rpt{k,1} = k;
    end
  else
    rpt{1} = 1:numel(tmpdata.trial);
    nrpt   = 1;
  end
  
  for rlop = 1:nrpt
    
    if dobvar % bvar
      for m = 1:ntap
        %---construct data-matrix
        for k = 1:size(cmbindx,1)
          dat = catnan(tmpdata.trial, cmbindx(k,:), rpt{rlop}, tap(m,:), nnans, dobvar);
          
          %---estimate autoregressive model
          switch cfg.method
            case 'biosig'
              [ar, rc, pe] = mvar(dat', cfg.order, cfg.mvarmethod);
              
              %---compute noise covariance
              tmpnoisecov     = pe(:,nchan*cfg.order+1:nchan*(cfg.order+1));
            case 'bsmart'
              [ar, tmpnoisecov] = armorf(dat, numel(rpt{rlop}), size(tmpdata.trial{1},2), cfg.order);
              ar = -ar; %convention is swapped sign with respect to biosig
              %FIXME check which is which: X(t) = A1*X(t-1) + ... + An*X(t-n) + E
              %the other is then X(t) + A1*X(t-1) + ... + An*X(t-n) = E
          end
          coeffs(rlop,:,k,:,j,m) = reshape(ar, [nchan*2 cfg.order]);
          
          %---rescale noisecov if necessary
          if dozscore, % FIX ME for bvar
            noisecov(rlop,:,k,:,j,m) = diag(datstd)*tmpnoisecov*diag(datstd);
          else
            noisecov(rlop,:,k,j,m) = reshape(tmpnoisecov,[1 4]);
          end
          dof(rlop,:,j) = numel(rpt{rlop});
        end
      end
    else % mvar
      for m = 1:ntap
        %---construct data-matrix
        dat = catnan(tmpdata.trial, chanindx, rpt{rlop}, tap(m,:), nnans, dobvar);
        
        %---estimate autoregressive model
        if dounivariate
          
          %---loop across the channels
          for p = 1:size(dat,1)
            
            switch cfg.method
              case 'biosig'
                [ar, rc, pe] = mvar(dat(p,:)', cfg.order, cfg.mvarmethod);
                
                %---compute noise covariance
                tmpnoisecov     = pe(:,cfg.order+1:(cfg.order+1));
              case 'bsmart'
                [ar, tmpnoisecov] = armorf(dat(p,:), numel(rpt{rlop}), size(tmpdata.trial{1},2), cfg.order);
                ar = -ar; %convention is swapped sign with respect to biosig
                %FIXME check which is which: X(t) = A1*X(t-1) + ... + An*X(t-n) + E
                %the other is then X(t) + A1*X(t-1) + ... + An*X(t-n) = E
            end
            coeffs(rlop,p,:,j,m) = reshape(ar, [1 cfg.order]);
            
            %---rescale noisecov if necessary
            if dozscore
              noisecov(rlop,p,j,m) = diag(datstd)*tmpnoisecov*diag(datstd);
            else
              noisecov(rlop,p,j,m) = tmpnoisecov;
            end
            dof(rlop,:,j) = numel(rpt{rlop});
          end
          
        else
          switch cfg.method
            case 'biosig'
              [ar, rc, pe] = mvar(dat', cfg.order, cfg.mvarmethod);
              
              %---compute noise covariance
              tmpnoisecov     = pe(:,nchan*cfg.order+1:nchan*(cfg.order+1));
            case 'bsmart'
              [ar, tmpnoisecov] = armorf(dat, numel(rpt{rlop}), size(tmpdata.trial{1},2), cfg.order);
              ar = -ar; %convention is swapped sign with respect to biosig
              %FIXME check which is which: X(t) = A1*X(t-1) + ... + An*X(t-n) + E
              %the other is then X(t) + A1*X(t-1) + ... + An*X(t-n) = E
          end
          coeffs(rlop,:,:,:,j,m) = reshape(ar, [nchan nchan cfg.order]);
          
          %---rescale noisecov if necessary
          if dozscore
            noisecov(rlop,:,:,j,m) = diag(datstd)*tmpnoisecov*diag(datstd);
          else
            noisecov(rlop,:,:,j,m) = tmpnoisecov;
          end
          dof(rlop,:,j) = numel(rpt{rlop});
        end %---dounivariate
        
      end %---tapers
    end
    
  end %---replicates
  
end %---tois
ft_progress('close');

%---create output-structure
mvardata          = [];

if ~dobvar && ~dounivariate && dojack
  mvardata.dimord = 'rptjck_chan_chan_lag';
elseif ~dobvar && ~dounivariate && keeprpt
  mvardata.dimord = 'rpt_chan_chan_lag';
elseif ~dobvar && ~dounivariate
  mvardata.dimord = 'chan_chan_lag';
  mvardata.label  = label;
  siz    = [size(coeffs) 1];
  coeffs = reshape(coeffs, siz(2:end));
  siz    = [size(noisecov) 1];
  if ~all(siz==1)
    noisecov = reshape(noisecov, siz(2:end));
  end
elseif dobvar
  mvardata.dimord = 'chancmb_lag';
  siz    = [size(coeffs) 1];
  coeffs = reshape(coeffs, [siz(2) * siz(3) siz(4) siz(5)]);
  siz    = [size(noisecov) 1];
  noisecov = reshape(noisecov, [siz(2) * siz(3) siz(4)]);
  mvardata.labelcmb = labelcmb;
elseif dounivariate
  mvardata.dimord = 'chan_lag';
  mvardata.label  = label;
  siz             = [size(coeffs) 1];
  coeffs          = reshape(coeffs, siz(2:end));
  siz    = [size(noisecov) 1];
  if ~all(siz==1)
    noisecov = reshape(noisecov, siz(2:end));
  end
end
mvardata.coeffs   = coeffs;
mvardata.noisecov = noisecov;
mvardata.dof      = dof;
if numel(cfg.toi)>1
  mvardata.time   = cfg.toi;
  mvardata.dimord = [mvardata.dimord,'_time'];
end
mvardata.fsampleorig = data.fsample;

switch cfg.output
  case 'parameters'
    % no output requested, do not re-compile time-series data
    
  case {'model' 'residual'}
    if keeprpt || dojack
      error('reconstruction of the residuals with keeprpt or dojack is not yet implemented');
    end
    
    dataout = keepfields(data, {'hdr','grad','fsample','trialinfo','label','cfg'});
    trial   = cell(1,numel(data.trial));
    time    = cell(1,numel(data.time));
    for k = 1:numel(data.trial)
      if strcmp(cfg.output, 'model')
        trial{k} = zeros(size(data.trial{k},1), size(data.trial{k},2)-cfg.order);
      else
        trial{k} = data.trial{k}(:, (cfg.order+1):end);
      end
      time{k}  = data.time{k}((cfg.order+1):end);
      for m = 1:cfg.order
        if dounivariate
          P = diag(mvardata.coeffs(:,m));
        else
          P = mvardata.coeffs(:,:,m);
        end
        
        if strcmp(cfg.output, 'residual')
          P = -P;
        end
        
        trial{k} = trial{k} + P * data.trial{k}(:,(cfg.order+1-m):(end-m));
      end
    end
    dataout.trial = trial;
    dataout.time  = time;
    
    cfg.coeffs   = mvardata.coeffs;
    cfg.noisecov = mvardata.noisecov;
    mvardata     = dataout; clear dataout;
  otherwise
    error('output ''%s'' is not supported', cfg.output);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance mvardata
ft_postamble history    mvardata
ft_postamble savevar    mvardata

%----------------------------------------------------
%subfunction to concatenate data with nans in between
function [datamatrix] = catnan(datacells, chanindx, trials, taper, nnans, dobvar)

nchan = length(chanindx);
nsmp  = cellfun('size', datacells, 2);
nrpt  = numel(trials);
sumsmp = cumsum([0 nsmp]);

%---initialize
datamatrix = nan(nchan, sum(nsmp) + nnans*(nrpt-1));

%---fill the matrix
for k = 1:nrpt
  if k==1
    begsmp = sumsmp(k) + 1;
    endsmp = sumsmp(k+1)  ;
  else
    begsmp = sumsmp(k)   + (k-1)*nnans + 1;
    endsmp = sumsmp(k+1) + (k-1)*nnans;
  end
  if ~dobvar && isempty(taper)
    datamatrix(:,begsmp:endsmp) = datacells{trials(k)}(chanindx,:);
  elseif ~dobvar && ~isempty(taper)
    % FIXME this will crash with variable data length and fixed length
    % taper
    datamatrix(:,begsmp:endsmp) = datacells{trials(k)}(chanindx,:).*taper(ones(nchan,1),:);
  elseif dobvar && isempty(taper)
    datamatrix(:,begsmp:endsmp) = datacells{trials(k)}(chanindx',:);
  elseif dobvar && ~isempty(taper)
    datamatrix(:,begsmp:endsmp) = datacells{trials(k)}(chanindx',:).*taper(ones(nchan,1),:);
  end
end

%------------------------------------------------------
%---subfunction to ensure that the first two input arguments are of double
% precision this prevents an instability (bug) in the computation of the
% tapers for MATLAB 6.5 and 7.0
function [tap] = double_dpss(a, b, varargin)
tap = dpss(double(a), double(b), varargin{:});
