function [cfg] = ft_neighbourplot(cfg, data)

% FT_NEIGHBOURPLOT visualizes neighbouring channels in a particular channel
% configuration. The positions of the channel are specified in a
% gradiometer or electrode configuration or from a layout.
%
% Use as
%   ft_neighbourplot(cfg)
% or as
%   ft_neighbourplot(cfg, data)
%
% where the configuration can contain
%   cfg.verbose       = 'yes' or 'no', if 'yes' then the plot callback will include text output
%   cfg.neighbours    = neighbourhood structure, see FT_PREPARE_NEIGHBOURS (optional)
% or one of the following options
%   cfg.layout        = filename of the layout, see FT_PREPARE_LAYOUT
%   cfg.elec          = structure with EEG electrode positions
%   cfg.grad          = structure with MEG gradiometer positions
%   cfg.elecfile      = filename containing EEG electrode positions
%   cfg.gradfile      = filename containing MEG gradiometer positions
%
% If cfg.neighbours is not defined or empty, this function will call
% FT_PREPARE_NEIGHBOURS to determine the channel neighbours. The
% following data fields may also be used by FT_PREPARE_NEIGHBOURS
%   data.elec     = structure with EEG electrode positions
%   data.grad     = structure with MEG gradiometer positions
%
% See also FT_PREPARE_NEIGHBOURS, FT_PREPARE_LAYOUT

% Copyright (C) 2011, J?rn M. Horschig, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

hasdata = nargin>1;
if hasdata, data = ft_checkdata(data); end

if isfield(cfg, 'neighbours')
  cfg.neighbours = cfg.neighbours;
elseif hasdata
  cfg.neighbours = ft_prepare_neighbours(cfg, data);
else
  cfg.neighbours = ft_prepare_neighbours(cfg);
end

if ~isfield(cfg, 'verbose')
  cfg.verbose = 'no';
elseif strcmp(cfg.verbose, 'yes')
  cfg.verbose = true;
end

% get the the grad or elec
if hasdata
  sens = ft_fetch_sens(cfg, data);
else
  sens = ft_fetch_sens(cfg);
end

% give some graphical feedback
if all(sens.chanpos(:,3)==0)
  % the sensor positions are already projected on a 2D plane
  proj = sens.chanpos(:,1:2);
else
  % use 3-dimensional data for plotting
  proj = sens.chanpos;
end
figure
axis equal
axis off
hold on;
for i=1:length(cfg.neighbours)
  this = cfg.neighbours(i);
  
  sel1 = match_str(sens.label, this.label);
  sel2 = match_str(sens.label, this.neighblabel);
  % account for missing sensors
  this.neighblabel = sens.label(sel2);
  for j=1:length(this.neighblabel)
    x1 = proj(sel1,1);
    y1 = proj(sel1,2);
    x2 = proj(sel2(j),1);
    y2 = proj(sel2(j),2);
    X = [x1 x2];
    Y = [y1 y2];
    if size(proj, 2) == 2
      line(X, Y, 'color', 'r');
    elseif size(proj, 2) == 3
      z1 = proj(sel1,3);
      z2 = proj(sel2(j),3);
      Z = [z1 z2];
      line(X, Y, Z, 'color', 'r');
    end
  end
end

% this is for putting the channels on top of the connections
for i=1:length(cfg.neighbours)
  this = cfg.neighbours(i);
  sel1 = match_str(sens.label, this.label);
  sel2 = match_str(sens.label, this.neighblabel);
  % account for missing sensors
  this.neighblabel = sens.label(sel2);
  if isempty(sel1)
    continue;
  end
  if size(proj, 2) == 2
    hs(i) = line(proj(sel1, 1), proj(sel1, 2),                                            ...
      'MarkerEdgeColor',  'k',                                        ...
      'MarkerFaceColor',  'k',                                        ...
      'Marker',           'o',                                        ...
      'MarkerSize',       .125*(2+numel(cfg.neighbours(i).neighblabel)^2), ...
      'UserData',         i,                                          ...
      'ButtonDownFcn',    @showLabelInTitle);
    
  elseif size(proj, 2) == 3
    hs(i) = line(proj(sel1, 1), proj(sel1, 2), proj(sel1, 3),                                ...
      'MarkerEdgeColor',  'k',                                        ...
      'MarkerFaceColor',  'k',                                        ...
      'Marker',           'o',                                        ...
      'MarkerSize',       .125*(2+numel(cfg.neighbours(i).neighblabel)^2), ...
      'UserData',         i,                        ...
      'ButtonDownFcn',    @showLabelInTitle);
  else
    error('Channel coordinates are too high dimensional');
  end
end
hold off;
title('[Click on a sensor to see its label]');

% store what is needed in UserData of figure
userdata.lastSensId = [];
userdata.cfg = cfg;
userdata.hs = hs;
set(gcf, 'UserData', userdata);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data

end % main function


function showLabelInTitle(gcbo, EventData, handles)

    userdata    = get(gcf, 'UserData');
    lastSensId  = userdata.lastSensId;
    cfg         = userdata.cfg;
    hs          = userdata.hs;
    curSensId   = get(gcbo, 'UserData');
    
    title(['Selected channel: ' cfg.neighbours(curSensId).label]);
    if cfg.verbose
      str = sprintf('%s, ', cfg.neighbours(curSensId).neighblabel{:});
      if length(str)>2
        % remove the last comma and space
        str = str(1:end-2);
      end
      fprintf('Selected channel %s, which has %d neighbours: %s\n', ...
        cfg.neighbours(curSensId).label, ...
        length(cfg.neighbours(curSensId).neighblabel), ...
        str);
    end
    
    set(hs(curSensId), 'MarkerFaceColor', 'g');
    set(hs(lastSensId), 'MarkerFaceColor', 'k');
    
    userdata.lastSensId = curSensId';    
    set(gcf, 'UserData', userdata);
  end
