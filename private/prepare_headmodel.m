function [headmodel, sens, cfg] = prepare_headmodel(cfg, data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that helps to prepare the electrodes/gradiometers and the
% volume conduction model. This is used in sourceanalysis and dipolefitting.
%
% This function will get the gradiometer/electrode definition and the volume
% conductor definition.
%
% Subsequently it will remove the gradiometers/electrodes that are not
% present in the data. Finally it with attach the gradiometers to a
% multi-sphere head model (if supplied) or attach the electrodes to
% the skin surface of a BEM head model.
%
% This function will return the electrodes/gradiometers in an order that is
% consistent with the order in cfg.channel, or in case that is empty in the
% order of the input electrode/gradiometer definition.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2004-2012, Robert Oostenveld
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

cfg = ft_checkconfig(cfg, 'forbidden', 'order');

% set the defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all');
cfg.siunits  = ft_getopt(cfg, 'siunits', 'no');   % yes/no, ensure that SI units are used consistently

hasdata = (nargin>1);

if hasdata
  % check if the input data is valid for this function
  data = ft_checkdata(data);
  % set the default for senstype depending on the data
  if isfield(data, 'grad')
    cfg.senstype = ft_getopt(cfg, 'senstype', 'meg');
  elseif isfield(data, 'elec')
    cfg.senstype = ft_getopt(cfg, 'senstype', 'eeg');
  elseif isfield(data, 'opto')
    cfg.senstype = ft_getopt(cfg, 'senstype', 'opto');
  else
    cfg.senstype = ft_getopt(cfg, 'senstype', []);
  end
end

% get the volume conduction model
if ischar(cfg.headmodel)
  headmodel = ft_read_headmodel(cfg.headmodel);
else
  % ensure that the volume conduction model is up-to-date
  headmodel = ft_datatype_headmodel(cfg.headmodel);
end

% get the gradiometer or electrode definition, these can be in the cfg or in the data
if hasdata
  sens = ft_fetch_sens(cfg, data);
else
  sens = ft_fetch_sens(cfg);
end

if istrue(cfg.siunits)
  % ensure that the geometrical units are in SI units
  sens       = ft_convert_units(sens,       'm', 'feedback', true);
  headmodel  = ft_convert_units(headmodel,  'm', 'feedback', true);
  if isfield(cfg, 'sourcemodel')
    cfg.sourcemodel = ft_convert_units(cfg.sourcemodel,  'm', 'feedback', true);
  end
else
  % ensure that the geometrical units are the same
  if isfield(cfg, 'sourcemodel') && isfield(cfg.sourcemodel, 'unit')
    % convert it to the units of the source model
    sens       = ft_convert_units(sens,       cfg.sourcemodel.unit, 'feedback', true);
    headmodel  = ft_convert_units(headmodel,  cfg.sourcemodel.unit, 'feedback', true);
  else
    % convert it to the units of the head model
    sens = ft_convert_units(sens, headmodel.unit, 'feedback', true);
  end
end

if hasdata && isfield(data, 'topolabel')
  % the data reflects a componentanalysis, where the topographic and the
  % timecourse labels are different
  cfg.channel = ft_channelselection(cfg.channel, data.topolabel);
elseif hasdata && isfield(data, 'label')
  % In the subsequent code, the matching channels in the sensor array and
  % in the configuration will be selected. To ensure that these channels
  % are also present in the data, update the configuration to match the data.
  cfg.channel = ft_channelselection(cfg.channel, data.label);
else
  % update the selected channels based on the electrode/gradiometer definition
  cfg.channel = ft_channelselection(cfg.channel, sens.label);
end

% ensure that these are a struct, which may be required in case configuration tracking is used
% FIXME this fails for combined EEG+MEG
% headmodel = struct(headmodel);
% sens      = struct(sens);

% the prepare_vol_sens function from the forwinv module does most of the actual work
[headmodel, sens] = ft_prepare_vol_sens(headmodel, sens, 'channel', cfg.channel);

% update the selected channels in the configuration
if iscell(sens)
  % this represents combined EEG, ECoG and/or MEG
  cfg.channel = {};
  for i=1:numel(sens)
    cfg.channel = cat(1, cfg.channel, sens{i}.label(:));
  end
else
  cfg.channel = sens.label;
end
