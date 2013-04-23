function [vol, sens, cfg] = prepare_headmodel(cfg, data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that helps to prepare the electrodes/gradiometers and the volume
% this is used in sourceanalysis and dipolefitting
%
% This function will get the gradiometer/electrode definition using
% FT_FETCH_SENS and the volume conductor definition using FT_FETCH_VOL
%
% Subsequently it will remove the gradiometers/electrodes that are not
% present in the data. Finally it with attach the gradiometers to a
% multi-sphere head model (if supplied) or attach the electrodes to
% a BEM head model.
%
% This function will return the electrodes/gradiometers in an order that is
% consistent with the order in cfg.channel, or in case that is empty in the
% order of the input electrode/gradiometer definition.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2004-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% set the defaults
if ~isfield(cfg, 'channel'),      cfg.channel = 'all';   end
if ~isfield(cfg, 'order'),        cfg.order = 10;        end % order of expansion for Nolte method; 10 should be enough for real applications; in simulations it makes sense to go higher
if ~isfield(cfg, 'sourceunits'),  cfg.sourceunits = 'cm';  end

if nargin<2
  data = [];
end

% get the volume conduction model
vol = ft_fetch_vol(cfg, data);

% get the gradiometer or electrode definition
sens = ft_fetch_sens(cfg, data);

% ensure that the units are the same
vol  = ft_convert_units(vol, cfg.sourceunits);
sens = ft_convert_units(sens, cfg.sourceunits);

if isfield(data, 'topolabel')
  % the data reflects a componentanalysis, where the topographic and the
  % timecourse labels are different
  cfg.channel = ft_channelselection(cfg.channel, data.topolabel);
elseif isfield(data, 'label')
  % In the subsequent code, the matching channels in the sensor array and
  % in the configuration will be selected. To ensure that these channels
  % are also present in the data, update the configuration to match the data.
  cfg.channel = ft_channelselection(cfg.channel, data.label);
else
  % update the selected channels based on the electrode/gradiometer definition
  cfg.channel = ft_channelselection(cfg.channel, sens.label);
end

% ensure that these are a struct, which may be required in case configuration tracking is used
vol  = struct(vol);
sens = struct(sens);

% the prepare_vol_sens function from the forwinv module does most of the actual work
[vol, sens] = ft_prepare_vol_sens(vol, sens, 'channel', cfg.channel, 'order', cfg.order);

% update the selected channels in the configuration
cfg.channel = sens.label;
