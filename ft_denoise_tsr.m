function data = ft_denoise_tsr(cfg, varargin)

% FT_DENOISE_TSR performs a regression analysis, using a (time-shifted set
% of) reference signal(s) as independent variable. It is a generic
% implementation of the method described by De Cheveigne ({REF}), or can be
% used to compute temporal-response-functions (see e.g. Lalor {REF}), or
% spatial filters based on canonical correlation (see Desain {REF})
%
% Use as
%   [dataout] = ft_denoise_tsr(cfg, data)
% or as
%   [dataout] = ft_denoise_tsr(cfg, data, refdata)
% where "data" is a raw data structure that was obtained with FT_PREPROCESSING. If
% you specify the additional input "refdata", the specified reference channels for
% the regression will be taken from this second data structure. This can be useful
% when reference-channel specific preprocessing needs to be done (e.g. low-pass
% filtering).
%
% The output structure dataout contains the denoised data in a format that is
% consistent with the output of FT_PREPROCESSING.
%
% The configuration should contain
%   cfg.refchannel = the channels used as reference signal (default = 'MEGREF')
%   cfg.channel    = the channels to be denoised (default = 'MEG')
%   cfg.method     = option specifying the criterion for the regression,
%                     'mvr', 'cca', 'pls', 'svd'
%   cfg.zscore     = standardise reference data prior to the regression (default = 'no')
%   cfg.updatesens = 'no' or 'yes' (default = 'yes')
%   cfg.perchannel
%   cfg.reflags
%
% if cfg.truncate is integer n > 1, n will be the number of singular values kept.
% if 0 < cfg.truncate < 1, the singular value spectrum will be thresholded at the
% fraction cfg.truncate of the explained variance.
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

% set the defaults
cfg.refchannel = ft_getopt(cfg, 'refchannel', 'MEGREF');
cfg.channel    = ft_getopt(cfg, 'channel',    'MEG');
cfg.truncate   = ft_getopt(cfg, 'truncate',   'no');
cfg.zscore     = ft_getopt(cfg, 'zscore',     'yes');
cfg.trials     = ft_getopt(cfg, 'trials',     'all', 1);
cfg.feedback   = ft_getopt(cfg, 'feedback',   'none');
cfg.updatesens = ft_getopt(cfg, 'updatesens', 'yes');
cfg.perchannel = ft_getopt(cfg, 'perchannel', 'yes');
cfg.reflags    = ft_getopt(cfg, 'reflags', 0);
cfg.method     = ft_getopt(cfg, 'method', 'mlr');
cfg.threshold  = ft_getopt(cfg, 'threshold', 0);
cfg.output     = ft_getopt(cfg, 'output', 'model');

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
[~, refdata] = rollback_provenance(cfg, refdata);  

% keep the requested channels from the data
tmpcfg  = keepfields(cfg, {'trials', 'showcallinfo' 'channel'});
data    = ft_selectdata(tmpcfg, varargin{1});
[cfg, data] = rollback_provenance(cfg, data);

% do the time shifting for the reference channel data
timestep = mean(diff(data.time{1}));
reflags = -round(cfg.reflags./timestep);
% the convention is to have a positive cfg.reflags defined as a delay of the ref w.r.t. the chan
% cellshift has an opposite convention with respect to the sign of the
% delay, hence the minus
if ~any(reflags==0)
  ft_error('the time lags for the reference data should at least include the sample 0');
end
fprintf('shifting the reference data\n');
refdata.trial = cellshift(refdata.trial, reflags, 2, [], 'overlap');

% center the data on lag 0
data.trial = cellshift(data.trial, 0, 2, [abs(min(reflags)) abs(max(reflags))], 'overlap');
data.time  = cellshift(data.time,  0, 2, [abs(min(reflags)) abs(max(reflags))], 'overlap');

% demean
fprintf('demeaning the data\n');
mu_refdata    = cellmean(refdata.trial, 2);
refdata.trial = cellvecadd(refdata.trial, -mu_refdata);
mu_data    = cellmean(data.trial, 2);
data.trial = cellvecadd(data.trial, -mu_data);

% zscore
if istrue(cfg.zscore)
  fprintf('zscoring the data\n');
  [refdata.trial, std_refdata] = cellzscore(refdata.trial, 2, 0);
  [data.trial,    std_data]    = cellzscore(data.trial, 2, 0);
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
    [Etmp, Dtmp] = multivariate_decomp(C(indx,indx), 1+(1:nref), 1, cfg.method, 1, cfg.threshold);
    beta_ref(k,:) = Etmp(2:end)./Etmp(1);
    rho(k)        = Dtmp;
  end
else
  ft_error('not yet implemented');
  %[E, D] = multivariate_decomp(C, 1:nchan, nchan+(1:nref), cfg.method, 1, cfg.threshold);  
end

dataout = keepfields(data, {'cfg' 'label' 'time' 'grad' 'elec' 'opto' 'trialinfo' 'fsample'});
switch cfg.output
  case 'model'
    dataout.trial = beta_ref*refdata.trial;
  case 'residual'
    dataout.trial = data.trial - beta_ref*refdata.trial;
end

if istrue(cfg.zscore)
  % unzscore the data
  dataout.trial = cellvecmult(dataout.trial, std_data);
end

weights.time = cfg.reflags;
weights.beta = beta_ref;
weights.rho  = rho;

dataout.weights = weights;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data

% rename the output variable to accomodate the savevar postamble
data = dataout;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
