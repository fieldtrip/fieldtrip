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
% Where the configuration can contain
%   cfg.verbose       = string, 'yes' or 'no', whether the function will print feedback text in the command window
%   cfg.neighbours    = neighbourhood structure, see FT_PREPARE_NEIGHBOURS (optional)
%   cfg.enableedit    = string, 'yes' or 'no', allows you to interactively add or remove edges between vertices (default = 'no')
%   cfg.visible       = string, 'on' or 'off' whether figure will be visible (default = 'on')
%   cfg.position      = location and size of the figure, specified as a vector of the form [left bottom width height]
%   cfg.renderer      = string, 'opengl', 'zbuffer', 'painters', see MATLAB Figure Properties. If this function crashes, you should try 'painters'.
%
% and either one of the following options
%   cfg.layout        = filename of the layout, see FT_PREPARE_LAYOUT
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad          = structure with gradiometer definition or filename, see FT_READ_SENS
%   cfg.opto          = structure with gradiometer definition or filename, see FT_READ_SENS
%
% If cfg.neighbours is not defined, this function will call
% FT_PREPARE_NEIGHBOURS to determine the channel neighbours. The
% following data fields may also be used by FT_PREPARE_NEIGHBOURS
%   data.elec         = structure with electrode positions
%   data.grad         = structure with gradiometer definition
%   data.opto         = structure with optode definition
%
% If cfg.neighbours is empty, no neighbouring sensors are assumed.
%
% Use cfg.enableedit to interactively add or remove edges in your own neighbour structure.
%
% See also FT_PREPARE_NEIGHBOURS, FT_PREPARE_LAYOUT

% Copyright (C) 2011, J?rn M. Horschig, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the data can be passed as input arguments or can be read from disk
hasdata = exist('data', 'var');

if hasdata
  % check if the input data is valid for this function
  data = ft_checkdata(data);
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed', {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed', {'optofile', 'opto'});
cfg = ft_checkconfig(cfg, 'renamed', {'newfigure', 'figure'});

% set the defaults
cfg.verbose    = ft_getopt(cfg, 'verbose', 'no');
cfg.enableedit = ft_getopt(cfg, 'enableedit', 'no');
cfg.visible    = ft_getopt(cfg, 'visible', 'on');
cfg.renderer   = ft_getopt(cfg, 'renderer', []); % let MATLAB decide on the default

if isfield(cfg, 'neighbours')
  cfg.neighbours = cfg.neighbours;
elseif hasdata
  cfg.neighbours = ft_prepare_neighbours(cfg, data);
else
  cfg.neighbours = ft_prepare_neighbours(cfg);
end

% get the the grad or elec
if hasdata
  sens = ft_fetch_sens(cfg, data);
else
  sens = ft_fetch_sens(cfg);
end

% insert sensors that are not in neighbourhood structure
if isempty(cfg.neighbours)
  nsel = 1:numel(sens.label);
else
  nsel = find(~ismember(sens.label, {cfg.neighbours.label}));
end

for i=1:numel(nsel)
  cfg.neighbours(end+1).label = sens.label{nsel(i)};
  cfg.neighbours(end).neighblabel = {};
end

[tmp, sel] = match_str(sens.label, {cfg.neighbours.label});
cfg.neighbours = cfg.neighbours(sel);

% give some graphical feedback
if all(sens.chanpos(:,3)==0)
  % the sensor positions are already projected on a 2D plane
  proj = sens.chanpos(:,1:2);
else
  % use 3-dimensional data for plotting
  proj = sens.chanpos;
end

% open a new figure with the specified settings
hf = open_figure(keepfields(cfg, {'figure', 'position', 'visible', 'renderer'}));

axis equal
axis vis3d
axis off
hold on

hl = [];
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
      hl(sel1, sel2(j)) = line(X, Y, 'color', 'r');
    elseif size(proj, 2) == 3
      z1 = proj(sel1,3);
      z2 = proj(sel2(j),3);
      Z = [z1 z2];
      hl(sel1, sel2(j)) = line(X, Y, Z, 'color', 'r');
    end
  end
end

% this is for putting the channels on top of the connections
hs = [];
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
      'MarkerSize',       .125*(2+numel(cfg.neighbours(i).neighblabel))^2, ...
      'UserData',         i,                                          ...
      'ButtonDownFcn',    @showLabelInTitle);

  elseif size(proj, 2) == 3
    hs(i) = line(proj(sel1, 1), proj(sel1, 2), proj(sel1, 3),                                ...
      'MarkerEdgeColor',  'k',                                        ...
      'MarkerFaceColor',  'k',                                        ...
      'Marker',           'o',                                        ...
      'MarkerSize',       .125*(2+numel(cfg.neighbours(i).neighblabel))^2, ...
      'UserData',         i,                        ...
      'ButtonDownFcn',    @showLabelInTitle);
  else
    ft_error('Channel coordinates are too high dimensional');
  end
end

hold off
title('[Click on a sensor to see its label]');

% store what is needed in UserData of figure
userdata.lastSensId = [];
userdata.cfg = cfg;
userdata.sens = sens;
userdata.hs = hs;
userdata.hl = hl;
userdata.quit = false;
hf = getparent(hf);
set(hf, 'UserData', userdata);

if istrue(cfg.enableedit)
  set(hf, 'CloseRequestFcn', @cleanup_cb);
  while ~userdata.quit
    uiwait(hf);
    userdata = get(hf, 'UserData');
  end
  cfg = userdata.cfg;

  hf = getparent(hf);
  delete(hf);
end

% remove SCALE and COMNT
desired = ft_channelselection({'all', '-SCALE', '-COMNT'}, {cfg.neighbours.label});

neighb_idx = ismember({cfg.neighbours.label}, desired);
cfg.neighbours = cfg.neighbours(neighb_idx);

% this is needed for the figure title
if isfield(cfg, 'dataname') && ~isempty(cfg.dataname)
  dataname = cfg.dataname;
elseif isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  dataname = cfg.inputfile;
elseif nargin>1
  dataname = arrayfun(@inputname, 2:nargin, 'UniformOutput', false);
else
  dataname = {};
end

% set the figure window title
if ~isempty(dataname)
  set(gcf, 'Name', sprintf('%d: %s: %s', double(gcf), mfilename, join_str(', ', dataname)));
else
  set(gcf, 'Name', sprintf('%d: %s', double(gcf), mfilename));
end
set(gcf, 'NumberTitle', 'off');

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous data
ft_postamble provenance
ft_postamble savefig

% add a menu to the figure, but only if the current figure does not have subplots
menu_fieldtrip(gcf, cfg, false);

if ~ft_nargout
  % don't return anything
  clear cfg
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showLabelInTitle(gcbo, EventData, handles)

userdata    = get(gcf, 'UserData');
lastSensId  = userdata.lastSensId;
cfg         = userdata.cfg;
hs          = userdata.hs;
curSensId   = get(gcbo, 'UserData');

if lastSensId == curSensId

  title('[Click on a sensor to see its label]');
  set(hs(curSensId), 'MarkerFaceColor', 'k');
  userdata.lastSensId = [];

elseif isempty(lastSensId) || ~istrue(cfg.enableedit)

  userdata.lastSensId = curSensId;
  if istrue(cfg.enableedit)
    title(['Selected channel: ' cfg.neighbours(curSensId).label ' click on another to (dis-)connect']);
  else
    title(['Selected channel: ' cfg.neighbours(curSensId).label]);
  end
  if istrue(cfg.verbose)
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

elseif istrue(cfg.enableedit)
  hl    = userdata.hl;
  sens  = userdata.sens;
  if all(sens.chanpos(:,3)==0)
    % the sensor positions are already projected on a 2D plane
    proj = sens.chanpos(:,1:2);
  else
    % use 3-dimensional data for plotting
    proj = sens.chanpos;
  end

  % find out whether they are connected
  connected1 = ismember(cfg.neighbours(curSensId).neighblabel, cfg.neighbours(lastSensId).label);
  connected2 = ismember(cfg.neighbours(lastSensId).neighblabel, cfg.neighbours(curSensId).label);

  if any(connected1) % then disconnect
    cfg.neighbours(curSensId).neighblabel(connected1) = [];
    cfg.neighbours(lastSensId).neighblabel(connected2) = [];
    title(['Disconnected channels ' cfg.neighbours(curSensId).label ' and ' cfg.neighbours(lastSensId).label]);
    delete(hl(curSensId, lastSensId));
    hl(curSensId, lastSensId) = 0;
    delete(hl(lastSensId, curSensId));
    hl(lastSensId, curSensId) = 0;
  else % then connect
    cfg.neighbours(curSensId).neighblabel{end+1} = cfg.neighbours(lastSensId).label;
    cfg.neighbours(lastSensId).neighblabel{end+1} = cfg.neighbours(curSensId).label;
    title(['Connected channels ' cfg.neighbours(curSensId).label ' and ' cfg.neighbours(lastSensId).label]);

    % draw new edge
    x1 = proj(curSensId,1);
    y1 = proj(curSensId,2);
    x2 = proj(lastSensId,1);
    y2 = proj(lastSensId,2);
    X = [x1 x2];
    Y = [y1 y2];
    if size(proj, 2) == 2
      hl(curSensId, lastSensId) = line(X, Y, 'color', 'r');
      hl(lastSensId, curSensId) = line(X, Y, 'color', 'r');
    elseif size(proj, 2) == 3
      z1 = proj(curSensId,3);
      z2 = proj(lastSensId,3);
      Z =[z1 z2];
      hl(curSensId, lastSensId) = line(X, Y, Z, 'color', 'r');
      hl(lastSensId, curSensId) = line(X, Y, Z, 'color', 'r');
    end

  end
  % draw nodes on top again
  delete(hs(curSensId));
  delete(hs(lastSensId));
  if size(proj, 2) == 2
    hs(curSensId) = line(proj(curSensId, 1), proj(curSensId, 2),                                            ...
      'MarkerEdgeColor',  'k',                                        ...
      'MarkerFaceColor',  'k',                                        ...
      'Marker',           'o',                                        ...
      'MarkerSize',       .125*(2+numel(cfg.neighbours(curSensId).neighblabel))^2, ...
      'UserData',         curSensId,                                          ...
      'ButtonDownFcn',    @showLabelInTitle);
    hs(lastSensId) = line(proj(lastSensId, 1), proj(lastSensId, 2),                                            ...
      'MarkerEdgeColor',  'k',                                        ...
      'MarkerFaceColor',  'k',                                        ...
      'Marker',           'o',                                        ...
      'MarkerSize',       .125*(2+numel(cfg.neighbours(lastSensId).neighblabel))^2, ...
      'UserData',         lastSensId,                                          ...
      'ButtonDownFcn',    @showLabelInTitle);

  elseif size(proj, 2) == 3
    hs(curSensId) = line(proj(curSensId, 1), proj(curSensId, 2), proj(curSensId, 3),                                ...
      'MarkerEdgeColor',  'k',                                        ...
      'MarkerFaceColor',  'k',                                        ...
      'Marker',           'o',                                        ...
      'MarkerSize',       .125*(2+numel(cfg.neighbours(curSensId).neighblabel))^2, ...
      'UserData',         curSensId,                        ...
      'ButtonDownFcn',    @showLabelInTitle);
    hs(lastSensId) = line(proj(lastSensId, 1), proj(lastSensId, 2), proj(lastSensId, 3),                                ...
      'MarkerEdgeColor',  'k',                                        ...
      'MarkerFaceColor',  'k',                                        ...
      'Marker',           'o',                                        ...
      'MarkerSize',       .125*(2+numel(cfg.neighbours(lastSensId).neighblabel))^2, ...
      'UserData',         lastSensId,                        ...
      'ButtonDownFcn',    @showLabelInTitle);
  else
    ft_error('Channel coordinates are too high dimensional');
  end

  if istrue(cfg.verbose)
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
  userdata.lastSensId = curSensId;
  userdata.hl = hl;
  userdata.hs = hs;
  userdata.cfg = cfg;
  set(gcf, 'UserData', userdata);
  return;
else
  % can never happen, so do nothing
end

set(gcf, 'UserData', userdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup_cb(h, eventdata)
userdata = get(h, 'UserData');
h   = getparent(h);
userdata.quit = true;
set(h, 'UserData', userdata);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end
