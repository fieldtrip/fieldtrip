function [selected] = ft_electrodeselection(cfg, elec)

% FT_ELECTRODESELECTION makes a selection of a predefined set of electrodes.
%
% Use as
%   selected = ft_electrodeselection(cfg, elec)
% where the input elec is a set of electrodes, for example from FT_READ_SENS or
% FT_ELECTRODEPLACEMENT, and cfg is a configuration structure that can contain the
% following options:
%
% If you want to make a selection on the basis of a specified list
%   cfg.channel   = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%
% If you want to make a graphical selection
%   cfg.headshape = a filename containing headshape, a structure containing a
%                   single triangulated boundary, or a Nx3 matrix with surface
%                   points
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_ELECTRODEPLACEMENT, FT_ELECTRODEREALIGN, FT_READ_SENS

% Copyright (C) 2026, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults                   % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init              % this will reset ft_warning and show the function help if nargin==0 and return an error
ft_preamble debug             % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar    elec % this reads the input data in case the user specified the cfg.inputfile option

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure the input data to be valid and update it to the latest standard
elec = ft_datatype_sens(elec);

% set the defaults
cfg.channel         = ft_getopt(cfg, 'channel', 'all');
cfg.headshape       = ft_getopt(cfg, 'headshape');
cfg.mesh            = ft_getopt(cfg, 'mesh');
cfg.axes            = ft_getopt(cfg, 'axes', 'yes');
cfg.elecstyle       = ft_getopt(cfg, 'elecstyle', {});
cfg.headshapestyle  = ft_getopt(cfg, 'headshapestyle', {});
cfg.meshstyle       = ft_getopt(cfg, 'meshstyle', {});

if ischar(cfg.headshape) && exist(cfg.headshape, 'file')
  ft_info('reading headshape from file %s\n', cfg.headshape);
  cfg.headshape = ft_read_headshape(cfg.headshape);
end

% set the default style using helper functions
cfg.elecstyle      = default_elecstyle(cfg.elecstyle, ~isempty(elec));
cfg.headshapestyle = default_headshapestyle(cfg.headshapestyle, ~isempty(cfg.headshape));
cfg.meshstyle      = default_meshstyle(cfg.meshstyle, ~isempty(cfg.mesh));

if ~isempty(cfg.headshape)
  % open a figure
  fig = open_figure(keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'}));
  set(fig, 'CloseRequestFcn',     @cb_quit);
  set(fig, 'WindowKeyPressFcn',   @cb_keyboard);
  set(fig, 'WindowButtonDownFcn', @cb_buttonpress);

  % make the initial selection
  cfg.channel = ft_channelselection(cfg.channel, elec.label);

  % create structure to be passed to gui
  opt.elec = elec;
  opt.headshape = cfg.headshape;
  opt.mesh = cfg.mesh;
  opt.showmesh = true;
  opt.elecstyle = cfg.elecstyle;
  opt.headshapestyle = cfg.headshapestyle;
  opt.meshstyle = cfg.meshstyle;
  opt.selected = ismember(elec.label, cfg.channel);
  opt.quit = false;
  setappdata(fig, 'opt', opt);

  cb_help(fig);
  while ~opt.quit
    cb_redraw(fig);
    uiwait(fig);
    opt = getappdata(fig, 'opt');
  end
  delete(fig);

  % make the final selection
  selected = make_selection(elec, opt.selected);
  % store the selection in the output cfg
  cfg.channel = elec.label(opt.selected);

else
  % only keep the desired channels, order them according to the users specification
  cfg.channel = ft_channelselection(cfg.channel, elec.label);
  [selchan, selsens] = match_str(cfg.channel, elec.label);

  % make the final selection
  selected = make_selection(elec, selsens);
end % if headshape

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble savevar selected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function elec = make_selection(elec, selected)
elec.label    = elec.label(selected);
elec.elecpos  = elec.elecpos(selected,:);
try, elec.elecori  = elec.elecori(selected,:);  end
try, elec.chantype = elec.chantype(selected,:); end
try, elec.chanunit = elec.chanunit(selected,:); end
try, elec.chanpos  = elec.chanpos (selected,:); end
try, elec.chanori  = elec.chanori (selected,:); end
try, elec.tra      = elec.tra(selected,:);      end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function style = default_elecstyle(style, present)
if present
  if iscell(style)
    % it should be a structure, not a cell-array
    style = ft_keyval2cfg(style);
  end
  % set the options for FT_PLOT_SENS
  style.label           = ft_getopt(style, 'label', 'off');
  style.chantype        = ft_getopt(style, 'chantype');
  style.orientation     = ft_getopt(style, 'orientation', false);
  % these have to do with the font
  style.fontcolor       = ft_getopt(style, 'fontcolor', 'k');  % default is black
  style.fontsize        = ft_getopt(style, 'fontsize',   get(0, 'defaulttextfontsize'));
  style.fontname        = ft_getopt(style, 'fontname',   get(0, 'defaulttextfontname'));
  style.fontweight      = ft_getopt(style, 'fontweight', get(0, 'defaulttextfontweight'));
  style.fontunits       = ft_getopt(style, 'fontunits',  get(0, 'defaulttextfontunits'));
  % this is for EEG electrodes
  style.elec            = ft_getopt(style, 'elec', false);
  style.elecshape       = ft_getopt(style, 'elecshape'); % default depends on the input, see below
  style.elecsize        = ft_getopt(style, 'elecsize');  % default depends on the input, see below
  style.headshape       = ft_getopt(style, 'headshape', []); % needed for elecshape 'disc'
else
  style = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function style = default_headshapestyle(style, present)
if present
  if iscell(style)
    % it should be a structure, not a cell-array
    style = ft_keyval2cfg(style);
  end
  % set the options for FT_PLOT_HEADSHAPE
  style.vertexcolor = ft_getopt(style, 'vertexcolor', 'none');
  style.facecolor   = ft_getopt(style, 'facecolor', 'skin');
  style.facealpha   = ft_getopt(style, 'facealpha', 1);
  style.edgecolor   = ft_getopt(style, 'edgecolor', 'none');
  style.vertexsize  = ft_getopt(style, 'vertexsize', 10);
  style.material    = ft_getopt(style, 'material', 'default');
  style.fidcolor    = ft_getopt(style, 'fidcolor', 'g');
  style.fidmarker   = ft_getopt(style, 'fidmarker', '*');
  style.fidlabel    = ft_getopt(style, 'fidlabel', true);
else
  style = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function style = default_meshstyle(style, present)
if present
  if iscell(style)
    % it should be a structure, not a cell-array
    style = ft_keyval2cfg(style);
  end
  % set the options for FT_PLOT_MESH
  style.vertexcolor   = ft_getopt(style, 'vertexcolor', 'none');
  style.facecolor     = ft_getopt(style, 'facecolor', [0.5 0.5 0.5]);
  style.edgecolor     = ft_getopt(style, 'edgecolor',   'none');
  style.faceindex     = ft_getopt(style, 'faceindex', false);
  style.vertexindex   = ft_getopt(style, 'vertexindex', false);
  style.vertexsize    = ft_getopt(style, 'vertexsize', 10);
  style.vertexmarker  = ft_getopt(style, 'vertexmarker', '.');
  style.facealpha     = ft_getopt(style, 'facealpha', 1);
  style.edgealpha     = ft_getopt(style, 'edgealpha', 1);
  style.edgelinewidth = ft_getopt(style, 'edgelinewidth', .5);
  style.material      = ft_getopt(style, 'material');
  style.surfaceonly   = ft_getopt(style, 'surfaceonly');
  style.clim          = ft_getopt(style, 'clim');
  style.alphalim      = ft_getopt(style, 'alphalim');
  style.alphamapping  = ft_getopt(style, 'alphamap', 'rampup');
  style.colormap      = ft_getopt(style, 'colormap');
  style.maskstyle     = ft_getopt(style, 'maskstyle', 'opacity');
  style.contour       = ft_getopt(style, 'contour', []);
  % these have to do with the font
  style.fontcolor       = ft_getopt(style, 'fontcolor', 'k');
  style.fontsize        = ft_getopt(style, 'fontsize',   get(0, 'defaulttextfontsize'));
  style.fontname        = ft_getopt(style, 'fontname',   get(0, 'defaulttextfontname'));
  style.fontweight      = ft_getopt(style, 'fontweight', get(0, 'defaulttextfontweight'));
  style.fontunits       = ft_getopt(style, 'fontunits',  get(0, 'defaulttextfontunits'));
  % these have to do with the contours
  style.contourcolor      = ft_getopt(style, 'contourcolor', 'k');
  style.contourlinewidth  = ft_getopt(style, 'contourlinewidth', 3);
  style.contourlinestyle  = ft_getopt(style, 'contourlinestyle', '-');
else
  style = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quit(h, eventdata)
fig = getparent(h);
opt = getappdata(h, 'opt');
opt.quit = true;
setappdata(fig, 'opt', opt);
uiresume;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)
h = getparent(h);
if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end
% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end
if isempty(key)
  % this happens if you press the apple key
  key = '';
end

% the following code is largely shared by FT_SOURCEPLOT, FT_VOLUMEREALIGN, FT_INTERACTIVEREALIGN, FT_MESHREALIGN, FT_ELECTRODEPLACEMENT
switch key
  case {'' 'shift+shift' 'alt-alt' 'control+control' 'command-0'}
    % do nothing
  case 'q'
    cb_quit(h);
  case 'h'
    cb_help(h);
  case 'm'
    opt = getappdata(h, 'opt');
    opt.showmesh = ~opt.showmesh;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
  otherwise
    % do nothing
end % switch key


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_help(h, eventdata)
disp('==================================================================================');
disp('0. Press "h" to show this help');
disp('1. Use the rotate option to change the viewpoint, selecting is disabled when rotating');
disp('2. Click on or close to an electrode to select or deselect it');
disp('3. Press "m" on the keyboard to toggle the visibility of the optional mesh')
disp('3. Press "q" on the keyboard or close the window when you are done');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonpress(h, eventdata)
h = getparent(h);
r = rotate3d(h);
if ~r.Enable
  cb_getposition(h);
end
cb_redraw(h);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_getposition(h, eventdata)
h   = getparent(h);
opt = getappdata(h, 'opt');

h2 = get(gca, 'children'); % get the object handles
iscorrect = false(size(h2));
for i=1:length(h2) % select the correct objects
  try
    pos = get(h2(i),'vertices');
    tri = get(h2(i),'faces');
    if ~isempty(opt.headshape) && isequal(opt.headshape.pos, pos) && isequal(opt.headshape.tri, tri)
      % it is the same object that the user has plotted before
      iscorrect(i) = true;
    elseif isempty(opt.headshape)
      % assume that it is the same object that the user has plotted before
      iscorrect(i) = true;
    end
  end
end
h2 = h2(iscorrect);
pos = select3d(h2)'; % enforce column direction
if ~isempty(pos)
  dist = opt.elec.elecpos;
  dist(:,1) = dist(:,1)-pos(1);
  dist(:,2) = dist(:,2)-pos(2);
  dist(:,3) = dist(:,3)-pos(3);
  dist = sqrt(sum(dist.^2,2));
  [mindist, indx] = min(dist);
  opt.selected(indx) = ~opt.selected(indx);
end
setappdata(h, 'opt', opt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)
h   = getparent(h);
opt = getappdata(h, 'opt');

% remember the current view
[az, el] = view();

% remember the current axis zoom
if ~isempty(get(gca, 'Children'))
  ax = axis();
else
  ax = [];
end

% redraw all geometrical objects
cla
axis vis3d
hold on

if ~isempty(opt.headshape)
  options = ft_cfg2keyval(opt.headshapestyle);
  ft_plot_headshape(opt.headshape, options{:}, 'axis', true);
end

if opt.showmesh && ~isempty(opt.mesh)
  options = ft_cfg2keyval(opt.meshstyle);
  if numel(opt.mesh)==length(opt.elec.label)
    % only plot the mesh for the selected electrodes
    ft_plot_mesh(opt.mesh(opt.selected), options{:});
  else
    ft_plot_mesh(opt.mesh, options{:});
  end
end


% plot the selected electrodes
ft_info('%d out of %d electrodes have been selected\n', sum(opt.selected), length(opt.selected));
options = ft_cfg2keyval(opt.elecstyle);
options = ft_setopt(options, 'facecolor', 'r');
ft_plot_sens(make_selection(opt.elec, opt.selected), options{:}); 

% plot the non-selected electrodes
options = ft_setopt(options, 'facecolor', [1 1 1]*0.8);
ft_plot_sens(make_selection(opt.elec, ~opt.selected), options{:});

% restore the current view
view(az, el);
if ~isempty(ax)
  axis(ax);
end

ft_headlight
rotate3d off % this is needed here, since ft_headlight switches it on
