function dataout = ft_denoise_tsr(cfg, varargin)

% FT_DENOISE_TSR performs a regression analysis, using a (time-shifted set
% of) reference signal(s) as independent variable. It is a generic
% implementation of the method described by De Cheveigne 
% (https://doi.org/10.1016/j.jneumeth.2007.06.003), or can be
% used to compute temporal-response-functions (see e.g. Crosse 
% (https://doi.org/10.3389/fnhum.2016.00604)), or
% spatial filters based on canonical correlation (see Thielen
% (https://doi.org/10.1371/journal.pone.0133797))
%
% Use as
%   [dataout] = ft_denoise_tsr(cfg, data)
% 
% or as
%   [dataout] = ft_denoise_tsr(cfg, data, refdata)
%
% where "data" is a raw data structure that was obtained with FT_PREPROCESSING. If
% you specify the additional input "refdata", the specified reference channels for
% the regression will be taken from this second data structure. This can be useful
% when reference-channel specific preprocessing needs to be done (e.g. low-pass
% filtering).
%
% The output structure dataout contains the denoised data in a format that is
% consistent with the output of FT_PREPROCESSING.
%
% The configuration options are:
%
%   cfg.refchannel         = the channels used as reference signal (default = 'MEGREF'), see FT_SELECTDATA
%   cfg.channel            = the channels to be denoised (default = 'all'), see FT_SELECTDATA 
%   cfg.method             = string, 'mlr', 'cca', 'pls', 'svd', option specifying the criterion for the regression
%                            (default = 'mlr')
%   cfg.reflags            = integer array, specifying temporal lags (in msec) by which to shift refchannel
%                            with respect to data channels
%   cfg.trials             = integer array, trials to be used in regression, see FT_SELECTDATA
%   cfg.testtrials         = cell array or string, trial indices to be used as test folds in a cross-validation scheme
%                            (numel(cfg.testrials == number of folds))
%   cfg.nfold              = scalar, indicating the number of test folds to
%                            use in a cross-validation scheme
%   cfg.standardiserefdata = string, 'yes' or 'no', whether or not to standardise reference data
%                            prior to the regression (default = 'no')
%   cfg.standardisedata    = string, 'yes' or 'no', whether or not to standardise dependent variable
%                            prior to the regression (default = 'no')
%   cfg.demeanrefdata      = string, 'yes' or 'no', whether or not to make
%                            reference data zero mean prior to the regression (default = 'no')
%   cfg.demeandata         = string, 'yes' or 'no', whether or not to make
%                            dependent variable zero mean prior to the regression (default = 'no')
%   cfg.threshold          = integer array, ([1 by 2] or [1 by numel(cfg.channel) + numel(cfg.reflags)]), 
%                            regularization or shrinkage ('lambda') parameter to be loaded on the diagonal of the
%                            penalty term (if cfg.method == 'mlrridge' or 'mlrqridge')
%   cfg.updatesens         = string, 'yes' or 'no' (default = 'yes')
%   cfg.perchannel         = string, 'yes' or 'no', or logical, whether or not to perform estimation of beta weights
%                            separately per channel
%   cfg.output             = string, 'model' or 'residual' (defaul = 'model'), 
%                            specifies what is outputed in .trial field in <dataout> 
%   cfg.performance        = string, 'Pearson' or 'r-squared' (default =
%                            'Pearson'), indicating what performance metric is outputed in .weights(k).performance
%                            field of <dataout> for the k-th fold     
%
% === cfg.threshold
% if cfg.threshold is 1 x 2 integer array, cfg.threshold(1) parameter scales uniformly
% in the dimension of predictor variable and cfg.threshold(2) in the space of
% response variable
%
%

% See also FT_PREPROCESSING, FT_DENOISE_SYNTHETIC, FT_DENOISE_PCA

% Copyright (c) 2008-2009, Jan-Mathijs Schoffelen, CCNi Glasgow
% Copyright (c) 2010-2011, Jan-Mathijs Schoffelen, DCCN Nijmegen
% Copyright (c) 2018, Jan-Mathijs Schoffelen
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

% UNDOCUMENTED OPTIONS (or possibly unused)
% cfg.testsamples 
% cfg.truncate   
% cfg.trials        
% 
% === cfg.truncate
% if cfg.truncate is integer n > 1, n will be the number of singular values kept.
% if 0 < cfg.truncate < 1, the singular value spectrum will be thresholded at the
% fraction cfg.truncate of the explained variance.

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance varargin
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'raw');
end

cfg.nfold       = ft_getopt(cfg, 'nfold',   1);
cfg.blocklength = ft_getopt(cfg, 'blocklength', 'trial');
cfg.testtrials  = ft_getopt(cfg, 'testtrials',  'all');
cfg.testsamples = ft_getopt(cfg, 'testsamples', 'all');
cfg.refchannel  = ft_getopt(cfg, 'refchannel', '');
cfg.reflags     = ft_getopt(cfg, 'reflags',    0); %this needs to be known for the folding

% set the rest of the defaults
cfg.channel            = ft_getopt(cfg, 'channel',            'all');
cfg.truncate           = ft_getopt(cfg, 'truncate',           'no');
cfg.standardiserefdata = ft_getopt(cfg, 'standardiserefdata', 'no');
cfg.standardisedata    = ft_getopt(cfg, 'standardisedata',    'no');
cfg.demeanrefdata      = ft_getopt(cfg, 'demeanrefdata',      'no');
cfg.demeandata         = ft_getopt(cfg, 'demeandata',         'no');
cfg.trials             = ft_getopt(cfg, 'trials',             'all', 1);
cfg.feedback           = ft_getopt(cfg, 'feedback',           'none');
cfg.updatesens         = ft_getopt(cfg, 'updatesens',         'yes');
cfg.perchannel         = ft_getopt(cfg, 'perchannel',         'yes');
cfg.method             = ft_getopt(cfg, 'method',             'mlr');
cfg.threshold          = ft_getopt(cfg, 'threshold',          0);
cfg.output             = ft_getopt(cfg, 'output',             'model');
cfg.performance        = ft_getopt(cfg, 'performance',        'Pearson');

if ~iscell(cfg.refchannel)
  cfg.refchannel = {cfg.refchannel};
end

if iscell(cfg.testtrials)
  % this has precedence above nfold
  cfg.nfold = numel(cfg.testtrials);
end

if cfg.nfold<=1
  dataout = ft_denoise_tsr_core(cfg, varargin{:});
else
  % do a cross validation
  if numel(varargin{1}.trial)>1 && ischar(cfg.blocklength) && isequal(cfg.blocklength, 'trial')
    if ~iscell(cfg.testtrials)
      % create sets of trial indices for the test data
      ntrl  = numel(varargin{1}.trial);
      edges = round(linspace(0,ntrl,cfg.nfold+1));
      indx  = randperm(ntrl);
      cfg.testtrials = cell(1,cfg.nfold);
      for k = 1:cfg.nfold
        cfg.testtrials{k} = indx((edges(k)+1):edges(k+1));
      end
    end
    
    testtrials = cfg.testtrials;
    tmp = cell(1,numel(testtrials));
    for k = 1:numel(testtrials)
      fprintf('estimating model for fold %d/%d\n', k, numel(testtrials));
      cfg.testtrials = testtrials{k};
      tmp{k}     = ft_denoise_tsr_core(cfg, varargin{:});
    end
    
    % create output data structure
    dataout = keepfields(tmp{1}, {'fsample' 'label'});
    for k = 1:numel(testtrials)
      tmp{k}.weights.trials = testtrials{k};
      
      dataout.trial(testtrials{k})   = tmp{k}.trial;
      dataout.time(testtrials{k})    = tmp{k}.time;
      dataout.weights(k)             = tmp{k}.weights;
      dataout.cfg.previous{k}        = tmp{k}.cfg;
      if isfield(tmp{k}, 'trialinfo')
        dataout.trialinfo(testtrials{k},:) = tmp{k}.trialinfo;
      end
    end
    
  elseif numel(varargin{1}.trial==1) ||(numel(varargin{1}.trial)>1 && ~ischar(cfg.blocklength))
    % concatenate into a single trial, with sufficient nan-spacing to
    % accommodate the shifting, and do a chunk-based folding
    error('not yet implemented');
  else
    error('incorrect specification of data and cfg.blocklength');
  end
  
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous varargin

ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

%-------------------------------------------------
function dataout = ft_denoise_tsr_core(cfg, varargin)

% create a separate structure for the reference data
tmpcfg  = keepfields(cfg, {'trials', 'showcallinfo'});
tmpcfg.channel = cfg.refchannel;
if numel(varargin)>1
  fprintf('selecting reference channel data from the second data input argument\n');
  refdata = ft_selectdata(tmpcfg, varargin{2});
else
  fprintf('selecting reference channel data from the first data input argument\n');
  refdata = ft_selectdata(tmpcfg, varargin{1});
end
[dum, refdata] = rollback_provenance(cfg, refdata);  

% keep the requested channels from the data
tmpcfg  = keepfields(cfg, {'trials', 'showcallinfo' 'channel'});
data    = ft_selectdata(tmpcfg, varargin{1});
[cfg, data] = rollback_provenance(cfg, data);

% deal with the specification of testtrials/testsamples, as per the
% instruction by the caller function, for cross-validation purposes
if ~ischar(cfg.testtrials) && ischar(cfg.testsamples) && isequal(cfg.testsamples, 'all')
  % subselect trials for testing
  usetestdata   = true;
  
  tmpcfg        = [];
  tmpcfg.trials = cfg.testtrials;
  testdata      = ft_selectdata(tmpcfg, data);
  testrefdata   = ft_selectdata(tmpcfg, refdata);
  tmpcfg.trials = setdiff(1:numel(data.trial), cfg.testtrials);
  data          = ft_selectdata(tmpcfg, data);
  refdata       = ft_selectdata(tmpcfg, refdata);
      
elseif ~ischar(cfg.testsamples) && ischar(cfg.testtrials) && isequal(cfg.testtrials, 'all')
  % subselect samples from a single trial for testing
  usetestdata = true;
elseif ischar(cfg.testtrials) && ischar(cfg.testsamples)
  % just a single fold, use all data for training and testing
  usetestdata = false;
else
  error('something wrong here');
end

% demean
if istrue(cfg.demeanrefdata)
  fprintf('demeaning the reference channels\n');
  mu_refdata    = cellmean(refdata.trial, 2);
  refdata.trial = cellvecadd(refdata.trial, -mu_refdata);
  if usetestdata
    mu_testrefdata    = cellmean(testrefdata.trial, 2);
    testrefdata.trial = cellvecadd(testrefdata.trial, -mu_testrefdata); 
  end
end
if istrue(cfg.demeandata)
  fprintf('demeaning the data channels\n');
  mu_data       = cellmean(data.trial, 2);
  data.trial    = cellvecadd(data.trial, -mu_data);
  if usetestdata
    mu_testdata    = cellmean(testdata.trial, 2);
    testdata.trial = cellvecadd(testdata.trial, -mu_testdata); 
  end
end

% standardise the data
if istrue(cfg.standardiserefdata)
  fprintf('standardising the reference channels \n');
  [refdata.trial, std_refdata] = cellzscore(refdata.trial, 2, 0);
end
if istrue(cfg.standardisedata)
  fprintf('standardising the data channels \n');
  [data.trial, std_data] = cellzscore(data.trial, 2, 0);
end

% do the time shifting for the reference channel data
ft_hastoolbox('cellfunction', 1);

timestep = mean(diff(data.time{1}));
reflags  = -round(cfg.reflags./timestep);
reflabel = refdata.label; % to be used later
% the convention is to have a positive cfg.reflags defined as a delay of the ref w.r.t. the chan
% cellshift has an opposite convention with respect to the sign of the
% delay, hence the minus
if ~any(reflags==0)
  ft_error('the time lags for the reference data should at least include the sample 0');
end
fprintf('shifting the reference data\n');
refdata.trial = cellshift(refdata.trial, reflags, 2, [], 'overlap');
refdata.time  = cellshift(data.time, 0, 2, [abs(min(reflags)) abs(max(reflags))], 'overlap');
refdata.label = repmat(refdata.label,numel(reflags),1);
for k = 1:numel(refdata.label)
  refdata.label{k} = sprintf('%s_shift%03d',refdata.label{k}, k);
end

% center the data on lag 0
data.trial = cellshift(data.trial, 0, 2, [abs(min(reflags)) abs(max(reflags))], 'overlap');
data.time  = cellshift(data.time,  0, 2, [abs(min(reflags)) abs(max(reflags))], 'overlap');

% only keep the trials that have > 0 samples
tmpcfg        = [];
tmpcfg.trials = find(cellfun('size',data.trial,2)>0);
data          = ft_selectdata(tmpcfg, data);
[cfg, data]   = rollback_provenance(cfg, data);
refdata       = ft_selectdata(tmpcfg, refdata);
[dum,refdata] = rollback_provenance(cfg, refdata);

% demean again, just to be sure
if istrue(cfg.demeanrefdata)
  fprintf('demeaning the reference channels\n');
  mu_refdata    = cellmean(refdata.trial, 2);
  refdata.trial = cellvecadd(refdata.trial, -mu_refdata);
end
if istrue(cfg.demeandata)
  fprintf('demeaning the data channels\n'); % the edges have been chopped off
  mu_data       = cellmean(data.trial, 2);
  data.trial    = cellvecadd(data.trial, -mu_data);
end

% compute the covariance
fprintf('computing the covariance\n');
nref  = size(refdata.trial{1},1);
nchan = numel(data.label);

C = nan(nchan,nchan);
C(1:nchan,1:nchan)        = nancov(data.trial,    data.trial, 1, 2, 1);
C(1:nchan,nchan+(1:nref)) = nancov(data.trial, refdata.trial, 1, 2, 1);
C(nchan+(1:nref),1:nchan) = C(1:nchan,nchan+(1:nref)).';
C(nchan+(1:nref),nchan+(1:nref)) = nancov(refdata.trial, refdata.trial, 1, 2, 1);

% compute the regression
if istrue(cfg.perchannel)
  beta_ref  = zeros(nchan, nref);
  rho = zeros(nchan,1);
  for k = 1:nchan
    indx = [k nchan+(1:nref)];
    [E, rho(k)]   = multivariate_decomp(C(indx,indx), 1+(1:nref), 1, cfg.method, 1, cfg.threshold);
    %beta_ref(k,:) = E(2:end)./E(1);
    beta_ref(k,:) = E(2:end);%./E(1);
  end
  %beta_ref = (diag(rho))*beta_ref; % scale with sqrt(rho), to get the proper scaling
  
else
  [E, rho]  = multivariate_decomp(C, 1:nchan, nchan+(1:nref), cfg.method, 1, cfg.threshold);  
  %beta_ref  = normc(E(nchan+(1:nref),:))';
  %beta_data = normc(E(1:nchan,:))';
  beta_ref  = E(nchan+(1:nref),:);
  beta_data = E(1:nchan,:);
end

% Unstandardise the data/refchannels and test data/refchannels
if istrue(cfg.standardiserefdata)
  std_refdata       = repmat(std_refdata, numel(cfg.reflags), 1);
  refdata.trial     = cellvecmult(refdata.trial, std_refdata);
  if exist('beta_data', 'var')
    beta_ref  = beta_ref'*diag(std_refdata);
  else
    beta_ref = diag(std_data)*beta_ref*diag(1./std_refdata);
  end
end

if istrue(cfg.standardisedata)
  data.trial     = cellvecmult(data.trial, std_data);
  if exist('beta_data', 'var')
    beta_data = diag(std_data)*beta_data;
  end
end

if usetestdata
  fprintf('shifting the reference data for the test data\n');
  testrefdata.trial = cellshift(testrefdata.trial, reflags, 2, [], 'overlap');
  testrefdata.time  = cellshift(testdata.time,  0, 2, [abs(min(reflags)) abs(max(reflags))], 'overlap');
  testrefdata.label = repmat(testrefdata.label,numel(reflags),1);
  for k = 1:numel(testrefdata.label)
    testrefdata.label{k} = sprintf('%s_shift%03d',testrefdata.label{k}, k);
  end
  
  % center the data on lag 0
  testdata.trial = cellshift(testdata.trial, 0, 2, [abs(min(reflags)) abs(max(reflags))], 'overlap');
  testdata.time  = cellshift(testdata.time,  0, 2, [abs(min(reflags)) abs(max(reflags))], 'overlap');
  
  % demean again, just to be sure
  if istrue(cfg.demeanrefdata)
    fprintf('demeaning the reference channels\n');
    mu_testrefdata    = cellmean(testrefdata.trial, 2);
    testrefdata.trial = cellvecadd(testrefdata.trial, -mu_testrefdata);
  end
  if istrue(cfg.demeandata)
    fprintf('demeaning the data channels\n'); % the edges have been chopped off
    mu_testdata       = cellmean(testdata.trial, 2);
    testdata.trial    = cellvecadd(testdata.trial, -mu_testdata);
  end

  % only keep the trials that have > 0 samples
  tmpcfg = [];
  tmpcfg.trials = find(cellfun('size',testdata.trial,2)>0);
  testdata     = ft_selectdata(tmpcfg, testdata);
  [dum,testdata] = rollback_provenance(cfg, testdata);
  testrefdata = ft_selectdata(tmpcfg, testrefdata);
  [dum,testrefdata] = rollback_provenance(cfg, testrefdata);

  predicted = beta_ref*testrefdata.trial;
  observed  = testdata.trial;
  time      = testdata.time;
else
  predicted = beta_ref*refdata.trial;
  observed  = data.trial;
  time      = data.time;
end

% create output data structure
dataout   = keepfields(data, {'cfg' 'label' 'grad' 'elec' 'opto' 'trialinfo' 'fsample'});
dataout.time = time;
switch cfg.output
  case 'model'
    dataout.trial = predicted;
  case 'residual'
    dataout.trial = observed - predicted;
end

% update the weights-structure
weights.time = cfg.reflags;
weights.rho  = rho;
if exist('beta_data', 'var')
  weights.unmixing = beta_data;
  weights.beta = beta_ref;
else
  % a per channel approach has been done, the beta weights reflect
  % (channelxtime-lag) -> reshape
  nref    = numel(cfg.refchannel);
  newbeta = zeros(size(beta_ref,1),size(beta_ref,2)./nref,nref);
  for k = 1:size(newbeta,3)
    newbeta(:,:,k) = beta_ref(:,k:nref:end);
  end
  weights.beta = newbeta;
  weights.reflabel = reflabel;
  weights.dimord   = 'chan_lag_refchan';
end

% Compute performance statistics
fprintf('Computing performance metric\n');
switch cfg.performance
  case 'Pearson'
    for k = 1:size(observed{1}, 1)
      tmp = nancov(cellcat(1, cellrowselect(observed,k), cellrowselect(predicted,k)), 1, 2, 1);
      weights.performance(k,1) = tmp(1,2)./sqrt(tmp(1,1).*tmp(2,2));
    end
  case 'r-squared' 
    tss = nansum(observed.^2, 2); % total sum of squares, testdata are already mean subtracted in l. 330
    rss = nansum((observed - predicted).^2, 2);  % sum of squared residual error
    % R-squared
    weights.performance = (tss-rss)./tss;
end

dataout.weights = weights;

