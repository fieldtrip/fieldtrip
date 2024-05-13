function [cfg] = ft_geometryplot(cfg)

% FT_GEOMETRYPLOT plots objects in 3D, such as sensors, headmodels, sourcemodels,
% headshapes, meshes, etcetera. It provides an easy-to-use wrapper for the
% corresponding FT_PLOT_XXX functions.
%
% Use as
%   ft_geometryplot(cfg)
% where the cfg structure contains the geometrical objects that have to be plotted
%   cfg.elec              = structure, see FT_READ_SENS
%   cfg.grad              = structure, see FT_READ_SENS
%   cfg.opto              = structure, see FT_READ_SENS
%   cfg.headshape         = structure, see FT_READ_HEADSHAPE
%   cfg.headmodel         = structure, see FT_PREPARE_HEADMODEL and FT_READ_HEADMODEL
%   cfg.sourcemodel       = structure, see FT_PREPARE_SOURCEMODEL
%   cfg.dipole            = structure, see FT_DIPOLEFITTING
%   cfg.mri               = structure, see FT_READ_MRI
%   cfg.mesh              = structure, see FT_PREPARE_MESH
%   cfg.axes              = string, 'yes' or 'no' (default = 'no')
%
% Furthermore, there are a number of general options
%   cfg.unit              = string, 'mm', 'cm', 'm' with the geometrical units (default depends on the data)
%   cfg.coordsys          = string, with the coordinate system (default depends on the data)
%   cfg.figure            = 'yes' or 'no', whether to open a new figure. You can also specify a figure handle from FIGURE, GCF or SUBPLOT. (default = 'yes')
%   cfg.figurename        = string, title of the figure window
%   cfg.position          = location and size of the figure, specified as [left bottom width height] (default is automatic)
%   cfg.renderer          = string, 'opengl', 'zbuffer', 'painters', see RENDERERINFO. The OpenGL renderer is required when using opacity (default = 'opengl')
%   cfg.title             = string, title of the plot
%
% You can specify the style with which the objects are displayed using
%   cfg.elecstyle         = cell-array or structure, see below
%   cfg.gradstyle         = cell-array or structure, see below
%   cfg.optostyle         = cell-array or structure, see below
%   cfg.headshapestyle    = cell-array or structure, see below
%   cfg.headmodelstyle    = cell-array or structure, see below
%   cfg.sourcemodelstyle  = cell-array or structure, see below
%   cfg.dipolestyle       = cell-array or structure, see below
%   cfg.mristyle          = cell-array or structure, see below
%   cfg.meshstyle         = cell-array or structure, see below
%
% For each of the xxxstyle options, you can specify a cell-array with key value pairs
% similar as in FT_INTERACTIVEREALIGN. These options will be passed on to the
% corresponding FT_PLOT_XXX function. You can also specify the options as a
% structure. For example, the following two specifications are equivalent
%   cfg.headshapestyle = {'facecolor', 'skin', 'edgecolor', 'none'};
% and
%   cfg.headshapestyle.facecolor = 'skin';
%   cfg.headshapestyle.edgecolor = 'none';
%
% In the figure with graphical user interface you will be able to adjust most of the
% settings that determine how the objects are displayed.
%
% See also PLOTTING, FT_SOURCEPLOT, FT_INTERACTIVEREALIGN

% Copyright (C) 2023, Robert Oostenveld
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
ft_preamble provenance

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% get the defaults
cfg.editfontsize  = ft_getopt(cfg, 'editfontsize',  12);
cfg.editfontunits = ft_getopt(cfg, 'editfontunits', 'points');     % inches, centimeters, normalized, points, pixels
cfg.renderer      = ft_getopt(cfg, 'renderer',      'opengl');
cfg.title         = ft_getopt(cfg, 'title',         []);
cfg.figurename    = ft_getopt(cfg, 'figurename',    []);
cfg.visible       = ft_getopt(cfg, 'visible',       'on');
cfg.coordsys      = ft_getopt(cfg, 'coordsys',      []);
cfg.unit          = ft_getopt(cfg, 'unit',          []);

% these are the geometrical objects that will be plotted
cfg.axes            = ft_getopt(cfg, 'axes', 'yes');
cfg.elec            = ft_getopt(cfg, 'elec', []);
cfg.grad            = ft_getopt(cfg, 'grad', []);
cfg.opto            = ft_getopt(cfg, 'opto', []);
cfg.headshape       = ft_getopt(cfg, 'headshape', []);
cfg.headmodel       = ft_getopt(cfg, 'headmodel', []);
cfg.sourcemodel     = ft_getopt(cfg, 'sourcemodel', []);
cfg.dipole          = ft_getopt(cfg, 'dipole', []);
cfg.mri             = ft_getopt(cfg, 'mri', []);
cfg.mesh            = ft_getopt(cfg, 'mesh', []);

% these are the options for plotting
cfg.elecstyle         = ft_getopt(cfg, 'elecstyle', {});
cfg.gradstyle         = ft_getopt(cfg, 'gradstyle', {});
cfg.optostyle         = ft_getopt(cfg, 'optostyle', {});
cfg.headshapestyle    = ft_getopt(cfg, 'headshapestyle', {});
cfg.headmodelstyle    = ft_getopt(cfg, 'headmodelstyle', {});
cfg.sourcemodelstyle  = ft_getopt(cfg, 'sourcemodelstyle', {});
cfg.dipolestyle       = ft_getopt(cfg, 'dipolestyle', {});
cfg.mristyle          = ft_getopt(cfg, 'mristyle', {});
cfg.meshstyle         = ft_getopt(cfg, 'meshstyle', {});

% set the defaults style using helper functions
cfg.elecstyle         = default_elecstyle(cfg.elecstyle,                ~isempty(cfg.elec));
cfg.gradstyle         = default_gradstyle(cfg.gradstyle,                ~isempty(cfg.grad));
cfg.optostyle         = default_optostyle(cfg.optostyle,                ~isempty(cfg.opto));
cfg.headshapestyle    = default_headshapestyle(cfg.headshapestyle,      ~isempty(cfg.headshape));
cfg.headmodelstyle    = default_headmodelstyle(cfg.headmodelstyle,      ~isempty(cfg.headmodel));
cfg.sourcemodelstyle  = default_sourcemodelstyle(cfg.sourcemodelstyle,  ~isempty(cfg.sourcemodel));
cfg.dipolestyle       = default_dipolestyle(cfg.dipolestyle,            ~isempty(cfg.dipole));
cfg.mristyle          = default_mristyle(cfg.mristyle,                  ~isempty(cfg.mri));
cfg.meshstyle         = default_meshstyle(cfg.meshstyle,                ~isempty(cfg.mesh));

if isempty(cfg.unit)
  % determine the common units
  tmp = {};
  for fn={'elec', 'grad', 'opto', 'headshape', 'headmodel', 'sourcemodel', 'dipole', 'mri', 'mesh'}
    if isfield(cfg, fn{1}) && isfield(cfg.(fn{1}), 'unit')
      tmp = [tmp {cfg.(fn{1}).unit}];
    end
  end
  if length(unique(tmp))==1
    cfg.unit = tmp{1};
  else
    ft_warning('inconsistent units in the input, converting all to ''mm''');
    cfg.unit = 'mm';
  end
end

if isempty(cfg.coordsys)
  % determine the common coordinate system
  tmp = {};
  for fn={'elec', 'grad', 'opto', 'headshape', 'headmodel', 'sourcemodel', 'dipole', 'mri', 'mesh'}
    if isfield(cfg, fn{1}) && isfield(cfg.(fn{1}), 'coordsys')
      tmp = [tmp {cfg.(fn{1}).coordsys}];
    end
  end
  if length(unique(tmp))==1
    cfg.coordsys = tmp{1};
  else
    ft_error('inconsistent coordinate systems in the input');
  end
end

% convert them all to the same units
for fn={'elec', 'grad', 'opto', 'headshape', 'headmodel', 'sourcemodel', 'dipole', 'mri', 'mesh'}
  if isfield(cfg, fn{1}) && ~isempty(cfg.(fn{1}))
    cfg.(fn{1}) = ft_convert_units(cfg.(fn{1}), cfg.unit);
  end
end

% open a new figure with the specified settings
fig = open_figure(keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'}));

setappdata(fig, 'cfg',         cfg); % this contains the data and most options
setappdata(fig, 'labels',      true);
setappdata(fig, 'grid',        true);
setappdata(fig, 'camlight',    true);
setappdata(fig, 'zoom',        0);    % between 0 and 1
setappdata(fig, 'elec',        true); % visible or not
setappdata(fig, 'grad',        true); % visible or not
setappdata(fig, 'opto',        true); % visible or not
setappdata(fig, 'headshape',   true); % visible or not
setappdata(fig, 'headmodel',   true); % visible or not
setappdata(fig, 'sourcemodel', true); % visible or not
setappdata(fig, 'dipole',      true); % visible or not
setappdata(fig, 'mri',         true); % visible or not
setappdata(fig, 'mesh',        true); % visible or not
setappdata(fig, 'axes',        true); % visible or not

% add the callbacks
set(fig, 'CloseRequestFcn',    @cb_quit);
set(fig, 'windowkeypressfcn',  @cb_keyboard);
set(fig, 'windowbuttondownfcn', @cb_buttonpress);

cb_creategui(fig);
cb_redraw(fig);
cb_help(cfg);
rotate3d on

% add a menu to the figure
menu_fieldtrip(fig, cfg, true);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble provenance

if ~ft_nargout
  % don't return anything
  clear cfg
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_creategui(h, eventdata, handles)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');

% move the axes so that they don't overlap with the GUI elements
set(gca, 'position', [0.15 0.30 0.50 0.60]);

% set the initial view
axis vis3d
view(3);

% navigation buttons on the bottom
uicontrol('tag', 'navx', 'parent', fig, 'style', 'pushbutton', 'string', '<<', 'callback', @cb_navigate);
uicontrol('tag', 'navx', 'parent', fig, 'style', 'pushbutton', 'string', '<',  'callback', @cb_navigate);
uicontrol('tag', 'navx', 'parent', fig, 'style', 'pushbutton', 'string', 'X',  'callback', @cb_navigate);
uicontrol('tag', 'navx', 'parent', fig, 'style', 'pushbutton', 'string', '>',  'callback', @cb_navigate);
uicontrol('tag', 'navx', 'parent', fig, 'style', 'pushbutton', 'string', '>>', 'callback', @cb_navigate);
uicontrol('tag', 'nav_', 'parent', fig, 'style', 'text', 'string', '', 'Visible', false);
uicontrol('tag', 'navy', 'parent', fig, 'style', 'pushbutton', 'string', '<<', 'callback', @cb_navigate);
uicontrol('tag', 'navy', 'parent', fig, 'style', 'pushbutton', 'string', '<',  'callback', @cb_navigate);
uicontrol('tag', 'navy', 'parent', fig, 'style', 'pushbutton', 'string', 'Y',  'callback', @cb_navigate);
uicontrol('tag', 'navy', 'parent', fig, 'style', 'pushbutton', 'string', '>',  'callback', @cb_navigate);
uicontrol('tag', 'navy', 'parent', fig, 'style', 'pushbutton', 'string', '>>', 'callback', @cb_navigate);
uicontrol('tag', 'nav_', 'parent', fig, 'style', 'text', 'string', '', 'Visible', false);
uicontrol('tag', 'navz', 'parent', fig, 'style', 'pushbutton', 'string', '<<', 'callback', @cb_navigate);
uicontrol('tag', 'navz', 'parent', fig, 'style', 'pushbutton', 'string', '<',  'callback', @cb_navigate);
uicontrol('tag', 'navz', 'parent', fig, 'style', 'pushbutton', 'string', 'Z',  'callback', @cb_navigate);
uicontrol('tag', 'navz', 'parent', fig, 'style', 'pushbutton', 'string', '>',  'callback', @cb_navigate);
uicontrol('tag', 'navz', 'parent', fig, 'style', 'pushbutton', 'string', '>>', 'callback', @cb_navigate);
% organize the buttons
ft_uilayout(fig, 'tag', 'nav.', 'units', 'normalized', 'width', 0.05, 'height', 0.05, 'hpos', 'auto', 'vpos', 0.05, 'backgroundcolor', [0.8 0.8 0.8], 'Visible', ~isempty(cfg.mri));

% control GUI elements on the right side
uicontrol('tag', 'elecbtn',        'parent', fig, 'style', 'checkbox',   'string', 'elec',                                     'value', getappdata(fig, 'elec'),        'callback', @cb_toggle, 'Visible', ~isempty(cfg.elec));
uicontrol('tag', 'gradbtn',        'parent', fig, 'style', 'checkbox',   'string', 'grad',                                     'value', getappdata(fig, 'grad'),        'callback', @cb_toggle, 'Visible', ~isempty(cfg.grad));
uicontrol('tag', 'optobtn',        'parent', fig, 'style', 'checkbox',   'string', 'opto',                                     'value', getappdata(fig, 'opto'),        'callback', @cb_toggle, 'Visible', ~isempty(cfg.opto));
uicontrol('tag', 'headshapebtn',   'parent', fig, 'style', 'checkbox',   'string', 'headshape',                                'value', getappdata(fig, 'headshape'),   'callback', @cb_toggle, 'Visible', ~isempty(cfg.headshape));
uicontrol('tag', 'headmodelbtn',   'parent', fig, 'style', 'checkbox',   'string', 'headmodel',                                'value', getappdata(fig, 'headmodel'),   'callback', @cb_toggle, 'Visible', ~isempty(cfg.headmodel));
uicontrol('tag', 'sourcemodelbtn', 'parent', fig, 'style', 'checkbox',   'string', 'sourcemodel',                              'value', getappdata(fig, 'sourcemodel'), 'callback', @cb_toggle, 'Visible', ~isempty(cfg.sourcemodel));
uicontrol('tag', 'dipolebtn',      'parent', fig, 'style', 'checkbox',   'string', 'dipole',                                   'value', getappdata(fig, 'dipole'),      'callback', @cb_toggle, 'Visible', ~isempty(cfg.dipole));
uicontrol('tag', 'mribtn',         'parent', fig, 'style', 'checkbox',   'string', 'mri',                                      'value', getappdata(fig, 'mri'),         'callback', @cb_toggle, 'Visible', ~isempty(cfg.mri));
uicontrol('tag', 'meshbtn',        'parent', fig, 'style', 'checkbox',   'string', 'mesh',                                     'value', getappdata(fig, 'mesh'),        'callback', @cb_toggle, 'Visible', ~isempty(cfg.mesh));
uicontrol('tag', 'axes1btn',       'parent', fig, 'style', 'checkbox',   'string', '3d axes',                                  'value', istrue(cfg.axes),               'callback', @cb_axes1);
uicontrol('tag', 'axes2btn',       'parent', fig, 'style', 'checkbox',   'string', 'figure axes',                              'value', getappdata(fig, 'axes'),        'callback', @cb_axes2);
uicontrol('tag', 'labelsbtn',      'parent', fig, 'style', 'checkbox',   'string', 'figure axes label',                        'value', getappdata(fig, 'labels'),      'callback', @cb_labels);
uicontrol('tag', 'gridbtn',        'parent', fig, 'style', 'checkbox',   'string', 'figure axes grid',                         'value', getappdata(fig, 'grid'),        'callback', @cb_grid);
uicontrol('tag', 'viewpointbtn',   'parent', fig, 'style', 'popup',      'string', 'default|top|bottom|left|right|front|back', 'value', 1,                              'callback', @cb_viewpoint);
uicontrol('tag', 'camlightbtn',    'parent', fig, 'style', 'popup',      'string', 'none|camlight|uniform',                    'value', 1,                              'callback', @cb_camlight);
uicontrol('tag', 'plotoptionsbtn', 'parent', fig, 'style', 'pushbutton', 'string', 'plot options ...',                                                                  'callback', @cb_plotoptions);

% remove the invisible GUI elements
delete(findall(fig, 'type', 'uicontrol', 'visible', false))

% organize the remaining ones
ft_uilayout(fig, 'tag', '.*btn', 'units', 'normalized', 'width', 0.25, 'height', 0.05, 'hpos', 0.75, 'vpos', 'auto', 'backgroundcolor', [0.8 0.8 0.8]);
ft_uilayout(fig, 'style', 'checkbox', 'backgroundcolor', get(fig, 'color'));

if isempty(cfg.coordsys) || strcmp(cfg.coordsys, 'unknown')
  ft_notice('the template coordinate system is unknown, selecting the viewpoint is not possible');
  ft_uilayout(fig, 'tag', 'viewpointbtn', 'Visible', 'off');
else
  [labelx, labely, labelz] = coordsys2label(cfg.coordsys, 2, 0);
  ft_notice('the template coordinate system is "%s"', cfg.coordsys);
  ft_info('the positive X-axis is pointing to %s', labelx);
  ft_info('the positive Y-axis is pointing to %s', labely);
  ft_info('the positive Z-axis is pointing to %s', labelz);
end
uiresume;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');

% disable all GUI elements
set(findall(fig, 'type', 'UIControl'), 'Enable', 'off')

% remember the current view
[az, el] = view();

% redraw all geometrical objects
cla
axis vis3d
hold on

% the options should be passed as a cell-array rather than as a structure
if getappdata(fig, 'elec') && ~isempty(cfg.elec)
  options = ft_cfg2keyval(cfg.elecstyle);
  ft_plot_sens(cfg.elec, options{:});
end

if getappdata(fig, 'grad') && ~isempty(cfg.grad)
  options = ft_cfg2keyval(cfg.gradstyle);
  ft_plot_sens(cfg.grad, options{:});
end

if getappdata(fig, 'opto') && ~isempty(cfg.opto)
  options = ft_cfg2keyval(cfg.optostyle);
  ft_plot_sens(cfg.opto, options{:});
end

if getappdata(fig, 'headshape') && ~isempty(cfg.headshape)
  options = ft_cfg2keyval(cfg.headshapestyle);
  ft_plot_headshape(cfg.headshape, options{:});
end

if getappdata(fig, 'headmodel') && ~isempty(cfg.headmodel)
  options = ft_cfg2keyval(cfg.headmodelstyle);
  ft_plot_headmodel(cfg.headmodel, options{:});
end

if getappdata(fig, 'sourcemodel') && ~isempty(cfg.sourcemodel)
  options = ft_cfg2keyval(cfg.sourcemodelstyle);
  ft_plot_mesh(cfg.sourcemodel.pos, options{:});
end

if getappdata(fig, 'dipole') && ~isempty(cfg.dipole)
  options = ft_cfg2keyval(cfg.dipolestyle);
  for i=1:size(cfg.dipole.pos,1)
    pos = cfg.dipole.pos(i,:);  % ndipoles * 3
    mom = cfg.dipole.mom(:,i);  % 3 * ndipoles
    ft_plot_dipole(pos, mom, options{:}, 'unit', cfg.unit);
  end
end

if getappdata(fig, 'mri') && ~isempty(cfg.mri)
  % this is required for a correct display
  cfg.mristyle = ft_setopt(cfg.mristyle, 'style', 'intersect');
  
  options = ft_cfg2keyval(cfg.mristyle);
  ft_plot_ortho(cfg.mri, options{:});

  % remove the labels from the MRI slices, FIXME there should be an option for this
  delete(findall(fig, 'Type', 'text', 'Tag', 'coordsyslabel_x'));
  delete(findall(fig, 'Type', 'text', 'Tag', 'coordsyslabel_y'));
  delete(findall(fig, 'Type', 'text', 'Tag', 'coordsyslabel_z'));
end

if getappdata(fig, 'mesh') && ~isempty(cfg.mesh)
  options = ft_cfg2keyval(cfg.meshstyle);
  ft_plot_mesh(cfg.mesh, options{:});
end

% restore the current view
view(az, el);

% update the figure based on the GUI elements
cb_camlight(h, []);
cb_axes1(h, []);
cb_axes2(h, []);
cb_labels(h, []);
cb_grid(h, []);

% re-enable all GUI elements
set(findall(fig, 'type', 'UIControl'), 'Enable', 'on')
uiresume;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_toggle(h, eventdata)
fig = getparent(h);
button = get(h, 'string');
if isvarname(button)
  % this switches the display of an object on or off
  setappdata(fig, button, get(h, 'value'));
end
cb_redraw(fig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_camlight(h, eventdata)
fig = getparent(h);
delete(findall(fig, 'Type', 'light'));  % remove all existing lights
switch get(findall(fig, 'tag', 'camlightbtn'), 'value')
  case 1 % none
    lighting none
  case 2 % camlight
    lighting gouraud
    camlight
  case 3 % uniform
    lighting gouraud
    uniformlight
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_axes1(h, eventdata)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');

if get(findall(fig, 'tag', 'axes1btn'), 'value')
  % FT_PLOT_AXES calls FT_PLOT_MESH, which does "axis off"
  ax = findall(fig, 'type', 'axes');
  prevaxis = ax.Visible;
  ft_plot_axes([], 'unit', cfg.unit, 'coordsys', cfg.coordsys, 'tag', 'axes1');
  % restore the figure axis in their previous state
  if prevaxis
    axis on
  else
    axis off
  end
else
  delete(findall(fig, 'tag', 'axes1'))
end
uiresume;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_axes2(h, eventdata)
fig = getparent(h);
if get(findall(fig, 'tag', 'axes2btn'), 'value')
  axis on
else
  axis off
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_labels(h, eventdata)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');
if get(findall(fig, 'tag', 'labelsbtn'), 'value')
  xlabel(sprintf('x (%s)', cfg.unit))
  ylabel(sprintf('y (%s)', cfg.unit))
  zlabel(sprintf('z (%s)', cfg.unit))
else
  xlabel('')
  ylabel('')
  zlabel('')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_grid(h, eventdata)
fig = getparent(h);
if get(findall(fig, 'tag', 'gridbtn'), 'value')
  grid on
else
  grid off
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_viewpoint(h, eventdata)
fig       = getparent(h);
coordsys  = getappdata(fig, 'coordsys');
% get the index of the option that was selected
val = get(h, 'value');
viewpoint = {'default', 'top', 'bottom', 'left', 'right', 'front', 'back'};
setviewpoint(gca, coordsys, viewpoint{val});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_help(h, eventdata)
fprintf('\n');
fprintf('==================================================================================\n');
fprintf('Press "h" to show this help.\n');
fprintf('Press "q" to quit.\n');
fprintf('Click and hold the left mouse button to rotate.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)
fig = getparent(h);

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(fig, 'userdata');
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

  case 'v' % camlight angle reset
    cb_camlight(h);

  otherwise
    % do nothing

end % switch key

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonpress(h, eventdata)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');

% this only works for the MRI slices, and only if they are orthogonal to the figure axes
obj = get(h, 'CurrentObject');
if isa(obj, 'matlab.graphics.primitive.Surface')
  ax = get(h, 'CurrentAxes');
  pos = get(ax, 'CurrentPoint');
  cfg.mristyle.location = mean(pos,1);
  setappdata(fig, 'cfg', cfg)
  cb_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_plotoptions(h, eventdata)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');

editfontsize  = cfg.editfontsize;
editfontunits = cfg.editfontunits;

code = [ ...
  sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') ...
  sprintf('%% Add or change options for on-the-fly preprocessing\n') ...
  sprintf('%% Use as cfg.xxxstyle.option = value                \n') ...
  sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') ...
  newline ...
  ];

for fn={'elecstyle', 'gradstyle', 'optostyle', 'headshapestyle', 'headmodelstyle', 'sourcemodelstyle', 'dipolestyle', 'mristyle', 'meshstyle'}
  style = cfg.(fn{1});
  if ~isempty(style)
    style = removelargefields(style);
    code = [code, printstruct(sprintf('cfg.%s', fn{1}), style), newline, newline];
  end
end

% make figure displaying the edit box
dialog = figure;
axis off

set(dialog, 'toolBar', 'none');
set(dialog, 'menuBar', 'none');
set(dialog, 'Position', get(fig, 'Position'));
set(dialog, 'Name', 'Edit plot options');
set(dialog, 'NumberTitle', 'off');

% add save button
uicontrol('parent', dialog, 'style', 'pushbutton', 'string', 'Cancel', 'callback', @cb_plotoptions_button);
uicontrol('parent', dialog, 'style', 'pushbutton', 'string', 'OK',     'callback', @cb_plotoptions_button);
ft_uilayout(dialog, 'style', 'pushbutton', 'units', 'normalized', 'hpos', 'auto', 'halign', 'right', 'vpos', 0.05)

% add edit box
editbox = uicontrol('parent', dialog, 'style', 'edit');
set(editbox, 'Units', 'normalized');
set(editbox, 'Position', [0.05 0.15 0.90 0.80]);
set(editbox, 'backgroundColor', [1 1 1]);
set(editbox, 'horizontalAlign', 'left');
set(editbox, 'max', 2);
set(editbox, 'min', 0);
set(editbox, 'FontName', 'Courier', 'FontSize', editfontsize, 'FontUnits', editfontunits);
set(editbox, 'string', code);

% add handle for the edit dialog to figure
setappdata(dialog, 'mainfigure', fig);
setappdata(dialog, 'editbox', editbox);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_plotoptions_button(h, eventdata)
dialogfigure = get(h, 'parent');
mainfigure = getappdata(dialogfigure, 'mainfigure');
button = get(h, 'string');

switch button
  case 'Cancel'
    close(dialogfigure);

  case 'OK'
    editbox = getappdata(dialogfigure, 'editbox');
    code = cellstr(get(editbox, 'string'));
    close(dialogfigure);

    % get rid of empty lines and white space
    skip = [];
    for iline = 1:numel(code)
      code{iline} = strtrim(code{iline});
      if isempty(code{iline})
        skip = [skip iline];
        continue
      end
      if code{iline}(1)=='%'
        skip = [skip iline];
        continue
      end
    end
    code(skip) = [];

    if ~isempty(code)
      iscorrect = startsWith(code, 'cfg.');
      if ~all(iscorrect)
        errordlg('Plot options must be specified as cfg.xxx = ...', 'Edit plot options', 'modal')
      end

      % eval the code line by line, this will update the local cfg variable
      cfg = getappdata(mainfigure, 'cfg');
      for icomm = 1:numel(code)
        eval([code{icomm} ';']);
      end
      setappdata(mainfigure, 'cfg', cfg);
      cb_redraw(mainfigure);
    end

end % switch button


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_navigate(h, eventdata)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');

% determine the size of the step
switch cfg.unit
  case 'mm'
    step = 1;
  case 'cm'
    step = 0.1;
  case 'm'
    step = 0.001;
end

switch get(h, 'string')
  case '<<'
    step = -10 * step;
  case '<'
    step = -1 * step;
  case '>'
    step = +1 * step;
  case '>>'
    step = +10 * step;
  case {'X', 'Y', 'Z'}
    % ask the user for the new location
    xstr = num2str(cfg.mristyle.location(1));
    ystr = num2str(cfg.mristyle.location(2));
    zstr = num2str(cfg.mristyle.location(3));
    answer = inputdlg({'x', 'y', 'z'}, 'location', [1 50], {xstr, ystr, zstr});
    if ~isempty(answer)
      cfg.mristyle.location(1) = str2double(answer{1});
      cfg.mristyle.location(2) = str2double(answer{2});
      cfg.mristyle.location(3) = str2double(answer{3});
      step = 0;
    end
end

dim = get(h, 'tag');
dim = dim(end); % x, y, or z

switch dim
  case 'x'
    cfg.mristyle.location(1) = cfg.mristyle.location(1) + step;
  case 'y'
    cfg.mristyle.location(2) = cfg.mristyle.location(2) + step;
  case 'z'
    cfg.mristyle.location(3) = cfg.mristyle.location(3) + step;
end

setappdata(fig, 'cfg', cfg);
cb_redraw(fig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quit(h, eventdata)
fig = getparent(h);
delete(fig)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = removelargefields(s)
fn = fieldnames(s);
for i=1:numel(fn)
  x = s.(fn{i});
  w = whos('x');
  b(i) = w.bytes;
end
skip = b>1000;
for i=1:numel(fn)
  x = s.(fn{i});
  if isnumeric(x) && all(size(x)>1)
    skip(i) = true;
  end
end
s = rmfield(s, fn(skip));

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
function style = default_gradstyle(style, present)
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
  % this is for MEG magnetometer and/or gradiometer arrays
  style.coil            = ft_getopt(style, 'coil', false);
  style.coilshape       = ft_getopt(style, 'coilshape'); % default depends on the input, see below
  style.coilsize        = ft_getopt(style, 'coilsize');  % default depends on the input, see below
else
  style = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function style = default_optostyle(style, present)
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
  style.fontcolor       = ft_getopt(style, 'fontcolor', 'k');
  style.fontsize        = ft_getopt(style, 'fontsize',   get(0, 'defaulttextfontsize'));
  style.fontname        = ft_getopt(style, 'fontname',   get(0, 'defaulttextfontname'));
  style.fontweight      = ft_getopt(style, 'fontweight', get(0, 'defaulttextfontweight'));
  style.fontunits       = ft_getopt(style, 'fontunits',  get(0, 'defaulttextfontunits'));
  % this is for NIRS optodes
  style.opto            = ft_getopt(style, 'opto', false);
  style.optoshape       = ft_getopt(style, 'optoshape');
  style.optosize        = ft_getopt(style, 'optosize');
  style.headshape       = ft_getopt(style, 'headshape'); % needed for optoshape 'disc'
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
function style = default_headmodelstyle(style, present)
if present
  if iscell(style)
    % it should be a structure, not a cell-array
    style = ft_keyval2cfg(style);
  end
  % set the options for FT_PLOT_HEADMODEL
  style.faceindex   = ft_getopt(style, 'faceindex', 'none');
  style.vertexindex = ft_getopt(style, 'vertexindex', 'none');
  style.vertexsize  = ft_getopt(style, 'vertexsize', 10);
  style.facecolor   = ft_getopt(style, 'facecolor', 'white');
  style.vertexcolor = ft_getopt(style, 'vertexcolor', 'none');
  style.edgecolor   = ft_getopt(style, 'edgecolor');
  style.facealpha   = ft_getopt(style, 'facealpha', 1);
  style.edgealpha   = ft_getopt(style, 'edgealpha', 1);
  style.surfaceonly = ft_getopt(style, 'surfaceonly');
else
  style = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function style = default_sourcemodelstyle(style, present)
if present
  if iscell(style)
    % it should be a structure, not a cell-array
    style = ft_keyval2cfg(style);
  end
  % set the options for FT_PLOT_MESH
  style.vertexcolor   = ft_getopt(style, 'vertexcolor', 'r');
  style.vertexindex   = ft_getopt(style, 'vertexindex', false);
  style.vertexsize    = ft_getopt(style, 'vertexsize', 10);
  style.vertexmarker  = ft_getopt(style, 'vertexmarker', '.');
else
  style = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function style = default_dipolestyle(style, present)
if present
  if iscell(style)
    % it should be a structure, not a cell-array
    style = ft_keyval2cfg(style);
  end
  % set the options for FT_PLOT_DIPOLE
  style.amplitudescale = ft_getopt(style, 'scale',     'none');
  style.color          = ft_getopt(style, 'color',     'r'); % can also be a RGB triplet
  style.alpha          = ft_getopt(style, 'alpha',      1);
  style.diameter       = ft_getopt(style, 'diameter',  'auto');
  style.length         = ft_getopt(style, 'length',    'auto');
  style.thickness      = ft_getopt(style, 'thickness', 'auto');
else
  style = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function style = default_mristyle(style, present)
if present
  if iscell(style)
    % it should be a structure, not a cell-array
    style = ft_keyval2cfg(style);
  end
  % set the options for FT_PLOT_ORTHO
  style.location            = ft_getopt(style, 'location', [0 0 0]);
  % style.orientation         = ft_getopt(style, 'orientation', eye(3));
  style.maskstyle           = ft_getopt(style, 'maskstyle', 'opacity');
  style.background          = ft_getopt(style, 'background');
  style.opacitylim          = ft_getopt(style, 'opacitylim');
  style.interpmethod        = ft_getopt(style, 'interpmethod', 'nearest');
  style.colormap            = ft_getopt(style, 'colormap');
  style.clim                = ft_getopt(style, 'clim');
  style.doscale             = ft_getopt(style, 'doscale', true);
  style.intersectcolor      = ft_getopt(style, 'intersectcolor', 'yrgbmyrgbm');
  style.intersectlinewidth  = ft_getopt(style, 'intersectlinewidth', 2);
  style.intersectlinestyle  = ft_getopt(style, 'intersectlinestyle');
  style.plotmarker          = ft_getopt(style, 'plotmarker');
  style.markersize          = ft_getopt(style, 'markersize', 'auto');
  style.markercolor         = ft_getopt(style, 'markercolor', 'w');
  style.facealpha           = ft_getopt(style, 'facealpha', 1);
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
  style.vertexcolor   = ft_getopt(style, 'vertexcolor', 'r');
  style.facecolor     = ft_getopt(style, 'facecolor', 'skin');
  style.edgecolor     = ft_getopt(style, 'edgecolor',   'k');
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
