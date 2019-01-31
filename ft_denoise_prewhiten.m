function [dataout] = ft_denoise_prewhiten(cfg, datain, noise)

% FT_DENOISE_PREWHITEN applies a spatial prewhitening operation to the data
% using the inverse noise covariance matrix. The consequence is that
% all channels are expressed in singnal-to-noise units, causing
% different channel types to be comparable. This ensures equal
% weighting in source estimation on data with different channel types.
%
% Use as
%   dataout = ft_denoise_prewhiten(cfg, datain, noise)
% where the datain is the original data from FT_PREPROCESSING, where
% noise should contain the estimated noise covariance from
% FT_TIMELOCKANALYSIS.
%
% The configuration structure can contain
%   cfg.channel     = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.split       = cell-array of channel types between which covariance is split, it can also be 'all' or 'no'
%   cfg.lambda      = scalar, or string
%   cfg.kappa       = scalar
%   cfg.tolerance   = scalar
%
% The channel selection relates to the channels that are pre-whitened
% using the same selection of channels in the noise covariance.
% All channels present in the input data structure will be present in the
% output, including trigger and
% other auxiliary channels.
%
% The lambda/kappa/tolerance values relate to how the inverse of the noise
% covariance is computed by FT_INV
%
% See also FT_DENOISE_SYNTHETIC, FT_DENOISE_PCA FT_INV

% Copyright (C) 2018-2019, Robert Oostenveld and Jan-Mathijs Schoffelen
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

ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    datain
ft_preamble provenance datain
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% get the defaults
cfg.channel = ft_getopt(cfg, 'channel', 'all');
cfg.split   = ft_getopt(cfg, 'split',   'all');
cfg.lambda  = ft_getopt(cfg, 'lambda',  0);
cfg.kappa   = ft_getopt(cfg, 'kappa',   []);
cfg.tol     = ft_getopt(cfg, 'tol',     []);
cfg.invmethod = ft_getopt(cfg, 'invmethod', 'tikhonov');

% ensure that the input data is correct, the next line is needed for a
% attempt correct detection of the data chanunit (with a hdr-field it fails
% for meggrad data)
if isfield(datain, 'hdr'), datain = rmfield(datain, 'hdr'); end

datain = ft_checkdata(datain, 'datatype', {'raw' 'timelock'}, 'haschantype', 'yes', 'haschanunit', 'yes'); 
noise  = ft_checkdata(noise,  'datatype',         'timelock', 'haschantype', 'yes', 'haschanunit', 'yes');

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'});
datain = ft_selectdata(tmpcfg, datain);
noise  = ft_selectdata(tmpcfg, noise);

% restore the provenance information
[cfg, datain] = rollback_provenance(cfg, datain);
[cfg, noise]  = rollback_provenance(cfg, noise);

if ~isfield(noise, 'cov')
  ft_error('noise covariance is not present');
end

% determine whether it is EEG and/or MEG data
haselec = isfield(datain, 'elec');
hasgrad = isfield(datain, 'grad');
% if haselec && hasgrad
%   ft_error('mixed covariance prewhitening is not supported for combined EEG/MEG');
% end

if isequal(cfg.split, 'no')
  chantype = {};
elseif isequal(cfg.split, 'all')
  chantype = unique(noise.chantype);
else
  chantype = cfg.split;
end

% zero out the off-diagonal elements for the specified channel types
for i=1:numel(chantype)
  sel = strcmp(noise.chantype, chantype{i});
  noise.cov(sel,~sel) = 0;
  noise.cov(~sel,sel) = 0;
end

% invert the noise covariance matrix
invnoise = ft_inv(noise.cov, 'lambda', cfg.lambda, 'kappa', cfg.kappa, 'tolerance', cfg.tol);
[U,S,V]  = svd(invnoise,'econ');

prewhiten             = [];
prewhiten.tra         = U*sqrt(S)*U';
prewhiten.labelold    = noise.label;
prewhiten.labelnew    = noise.label;
prewhiten.chantypeold = noise.chantype;
prewhiten.chantypenew = noise.chantype;
prewhiten.chanunitold = noise.chanunit;
prewhiten.chanunitnew = repmat({'snr'}, size(noise.chantype));

% apply the projection to the data
dataout = ft_apply_montage(datain, prewhiten, 'keepunused', 'yes');

if hasgrad
  % the gradiometer structure needs to be updated to ensure that the forward model remains consistent with the data
  dataout.grad = ft_apply_montage(datain.grad, prewhiten, 'balancename', 'prewhiten');
end

if haselec
  % the electrode structure needs to be updated to ensure that the forward model remains consistent
  dataout.elec = ft_apply_montage(datain.elec, prewhiten, 'balancename', 'prewhiten');
end

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

