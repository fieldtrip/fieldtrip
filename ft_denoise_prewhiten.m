function [dataout] = ft_denoise_prewhiten(cfg, datain, noise)

% FT_DENOISE_PREWHITEN applies a prewhitening operation to the data
% using the inverse noise covariance matrix. The consequence is that
% all channels are expressed in singnal-to-noise units, causing
% different channel types to be comparable. This ensures equal
% weithing in source estimation on data with different channel types.
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
%
% See also FT_DENOISE_SYNTHETIC, FT_DENOISE_PCA

% Copyright (C) 2018, Robert Oostenveld
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
cfg.mixed   = ft_getopt(cfg, 'channel', 'all');

% ensure that the input data is correct
datain = ft_checkdata(datain, 'datatype', 'raw', 'isnirs', 'no');
noise  = ft_checkdata(noise, 'datatype', 'timelock');

if ~isfield(noise, 'cov')
  ft_error('noise covariance is not present');
end

% determine whether it is EEG and/or MEG data
haselec = isfield(datain, 'elec');
hasgrad = isfield(datain, 'grad');
if haselec && hasgrad && istrue(cfg.mixed)
  ft_error('mixed covariance prewhitening is not supported for combined EEG/MEG');
end

% get the intersection
cfg.channel = ft_channelselection(cfg.channel) datain.label);
cfg.channel = ft_channelselection(cfg.channel) noise.label);

% select the overlapping channels
tmpcfg = keepfields(cfg, 'channel');
datain = ft_selectdata(tmpcfg, datain);
[~, datain] = rollback_provenance(cfg, datain);
noise = ft_selectdata(tmpcfg, noise);
[~, noise] = rollback_provenance(cfg, noise);

% do a sanity check on the channel types
chantype = unique(datain.chantype);
for i=1:numel(chantype)
  sel = strcmp(datain.chantype, chantype{i});
  tmp = datain.chanunit(sel);
  if ~all(strcmp(tmp, tmp{1});
    ft_warning('not all %s channels have the same units', chantype{i});
  end
end

% make the channels of the noise covariance consistent with the data
[datsel, noisesel] = match_str(datain.label, noise.label);
noisecov = noise.cov(noisesel,noisesel);

if isequal(cfg.split, 'no')
  chantype = {};
elseif isequal(cfg.split, 'all')
  chantype = unique(datain.chantype);
else
  chantype = cfg.split;
end

% zero out the off-diagonal elements for the specified channel types
for i=1:numel(chantype)
  sel = strcmp(datain.chantype, chantype{i});
  noisecov(sel,~sel) = 0;
  noisecov(~sel,sel) = 0;
end

prewhiten = [];
prewhiten.tra = inv(noisecov); % FIXME add regularization
prewhiten.labelold = datain.label;
prewhiten.labelnew = datain.label;
prewhiten.chantypeold = datain.chantype;
prewhiten.chantypenew = datain.chantype;
prewhiten.chanunitold = datain.chanunit;
prewhiten.chanunitnew = repmat({'snr', size(datain.chantype));

% apply the projection to the data
dataout = ft_apply_montage(datain, prewhiten);

if isgrad
  % the gradiometer structure needs to be updated to ensure that the forward model remains consistent with the data
  dataout.grad = ft_apply_montage(datain.grad, prewhiten, 'balancename', 'prewhiten');
end

if iselec
  % the electrode structure needs to be updated to ensure that the forward model remains consistent
  dataout.elec = ft_apply_montage(datain.elec, prewhiten, 'balancename', 'prewhiten');
end

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   datain
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

