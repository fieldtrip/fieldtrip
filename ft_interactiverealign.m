function [cfg] = ft_interactiverealign(cfg)

% FT_INTERACTIVEREALIGN allows the user to interactively translate, rotate and scale an
% individual geometrical object to a template geometrical object. It can for example be used
% to align EEG electrodes to a model of the scalp surface.
%
% Use as
%   [cfg] = ft_interactiverealign(cfg)
%
% The configuration structure should contain the individuals geometrical object that
% has to be realigned
%   cfg.individual.elec           = structure, see FT_READ_SENS
%   cfg.individual.grad           = structure, see FT_READ_SENS
%   cfg.individual.opto           = structure, see FT_READ_SENS
%   cfg.individual.headmodel      = structure, see FT_PREPARE_HEADMODEL
%   cfg.individual.headshape      = structure, see FT_READ_HEADSHAPE
%   cfg.individual.mri            = structure, see FT_READ_MRI
%   cfg.individual.mesh           = structure, see FT_PREPARE_MESH
% You can specify the style with which the objects are displayed using
%   cfg.individual.headmodelstyle = 'vertex', 'edge', 'surface' or 'both' (default = 'edge')
%   cfg.individual.headshapestyle = 'vertex', 'edge', 'surface' or 'both' (default = 'vertex')
%
% The configuration structure should also contain the geometrical object of a
% template that serves as target
%   cfg.template.axes             = string, 'yes' or 'no' (default = 'no')
%   cfg.template.elec             = structure, see FT_READ_SENS
%   cfg.template.grad             = structure, see FT_READ_SENS
%   cfg.template.opto             = structure, see FT_READ_SENS
%   cfg.template.headmodel        = structure, see FT_PREPARE_HEADMODEL
%   cfg.template.headshape        = structure, see FT_READ_HEADSHAPE
%   cfg.template.mri              = structure, see FT_READ_MRI
%   cfg.template.mesh             = structure, see FT_PREPARE_MESH
% You can specify the style with which the objects are displayed using
%   cfg.template.headmodelstyle   = 'vertex', 'edge', 'surface' or 'both' (default = 'edge')
%   cfg.template.headshapestyle   = 'vertex', 'edge', 'surface' or 'both' (default = 'vertex')
%
% You can specify one or multiple individual objects which will all be realigned and
% one or multiple template objects.
%
% See also FT_VOLUMEREALIGN, FT_ELECTRODEREALIGN, FT_DETERMINE_COORDSYS,
% FT_READ_SENS, FT_READ_HEADMODEL, FT_READ_HEADSHAPE

% Undocumented options, primarily to support FT_DEFACEVOLUME
%   cfg.showalpha
%   cfg.showlight
%   cfg.showmaterial
%   cfg.showapply
%   cfg.rotate
%   cfg.scale
%   cfg.translate
%   cfg.transformorder

% Copyright (C) 2008, Vladimir Litvak
% Copyright (C) 2022-2024, Robert Oostenveld
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
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', {'individual', 'template'});
cfg.individual = ft_checkconfig(cfg.individual, 'renamed', {'vol', 'headmodel'});
cfg.individual = ft_checkconfig(cfg.individual, 'renamed', {'volstyle', 'headmodelstyle'});
cfg.template   = ft_checkconfig(cfg.template, 'renamed', {'vol', 'headmodel'});
cfg.template   = ft_checkconfig(cfg.template, 'renamed', {'volstyle', 'headmodelstyle'});

% allow to pass an initial transformation, this is used by FT_DEFACEVOLUME
cfg.rotate         = ft_getopt(cfg, 'rotate', [0 0 0]);
cfg.scale          = ft_getopt(cfg, 'scale', [1 1 1]);
cfg.translate      = ft_getopt(cfg, 'translate', [0 0 0]);
cfg.transformorder = ft_getopt(cfg, 'transformorder', {'rotate', 'translate', 'scale'});

% get the other options
cfg.showalpha                 = ft_getopt(cfg, 'showalpha'); % default is set below
cfg.showlight                 = ft_getopt(cfg, 'showlight', 'yes');
cfg.showmaterial              = ft_getopt(cfg, 'showlight', 'yes');
cfg.showapply                 = ft_getopt(cfg, 'showapply', 'yes');
cfg.unit                      = ft_getopt(cfg, 'unit', 'mm');
cfg.individual.elec           = ft_getopt(cfg.individual, 'elec', []);
cfg.individual.elecstyle      = ft_getopt(cfg.individual, 'elecstyle', {}); % key-value pairs
cfg.individual.grad           = ft_getopt(cfg.individual, 'grad', []);
cfg.individual.gradstyle      = ft_getopt(cfg.individual, 'gradstyle', {}); % key-value pairs
cfg.individual.opto           = ft_getopt(cfg.individual, 'opto', []);
cfg.individual.optostyle      = ft_getopt(cfg.individual, 'optostyle', {}); % key-value pairs
cfg.individual.headshape      = ft_getopt(cfg.individual, 'headshape', []);
if ~isempty(cfg.individual.headshape) && isfield(cfg.individual.headshape, 'tri')
  cfg.individual.headshapestyle = ft_getopt(cfg.individual, 'headshapestyle', 'surface');
else
  cfg.individual.headshapestyle = ft_getopt(cfg.individual, 'headshapestyle', 'vertex');
end
cfg.individual.headmodel      = ft_getopt(cfg.individual, 'headmodel', []);
cfg.individual.headmodelstyle = ft_getopt(cfg.individual, 'headmodelstyle', 'edge');
cfg.individual.mri            = ft_getopt(cfg.individual, 'mri', []);
cfg.individual.mristyle       = ft_getopt(cfg.individual, 'mristyle', {});
cfg.individual.mesh           = ft_getopt(cfg.individual, 'mesh', []);
cfg.individual.meshstyle      = ft_getopt(cfg.individual, 'meshstyle', {});

cfg.template.axes             = ft_getopt(cfg.template, 'axes', 'no');
cfg.template.elec             = ft_getopt(cfg.template, 'elec', []);
cfg.template.elecstyle        = ft_getopt(cfg.template, 'elecstyle', {}); % key-value pairs
cfg.template.grad             = ft_getopt(cfg.template, 'grad', []);
cfg.template.gradstyle        = ft_getopt(cfg.template, 'gradstyle', {}); % key-value pairs
cfg.template.opto             = ft_getopt(cfg.template, 'opto', []);
cfg.template.optostyle        = ft_getopt(cfg.template, 'optostyle', {}); % key-value pairs
cfg.template.headshape        = ft_getopt(cfg.template, 'headshape', []);
if ~isempty(cfg.template.headshape) && isfield(cfg.template.headshape, 'tri')
  cfg.template.headshapestyle   = ft_getopt(cfg.template, 'headshapestyle', 'surface');
else
  cfg.template.headshapestyle   = ft_getopt(cfg.template, 'headshapestyle', 'vertex');
end
cfg.template.headmodel        = ft_getopt(cfg.template, 'headmodel', []);
cfg.template.headmodelstyle   = ft_getopt(cfg.template, 'headmodelstyle', 'edge');
cfg.template.mri              = ft_getopt(cfg.template, 'mri', []);
cfg.template.mristyle         = ft_getopt(cfg.template, 'mristyle', {});
cfg.template.mesh             = ft_getopt(cfg.template, 'mesh', []);
cfg.template.meshstyle        = ft_getopt(cfg.template, 'meshstyle', {});

% convert the string that describes the style to a cell-array
cfg.template.headshapestyle   = updatestyle(cfg.template.headshapestyle);
cfg.individual.headshapestyle = updatestyle(cfg.individual.headshapestyle);
cfg.template.headmodelstyle   = updatestyle(cfg.template.headmodelstyle);
cfg.individual.headmodelstyle = updatestyle(cfg.individual.headmodelstyle);

% ensure that these are keyval cell-arrays, not cfg structures
for fn1={'individual', 'template'}
  fn1 = fn1{1};
  for fn2={'elecstyle', 'gradstyle', 'headshapestyle', 'headmodelstyle', 'mristyle', 'meshstyle'}
    fn2 = fn2{1};
    style = cfg.(fn1).(fn2);
    if isstruct(style)
      cfg.(fn1).(fn2) = ft_cfg2keyval(style);
    end
  end
end

if ~isempty(cfg.individual.headshape) && isfield(cfg.individual.headshape, 'color')
  if isfield(cfg.individual.headshape, 'tri') && size(cfg.individual.headshape.tri,1)==size(cfg.individual.headshape.color,1)
    cfg.individual.headshapestyle = ft_setopt(cfg.individual.headshapestyle, 'facecolor', cfg.individual.headshape.color);
  elseif size(cfg.individual.headshape.pos,1)==size(cfg.individual.headshape.color,1)
    cfg.individual.headshapestyle = ft_setopt(cfg.individual.headshapestyle, 'vertexcolor', cfg.individual.headshape.color);
  end
end

if ~isempty(cfg.template.headshape) && isfield(cfg.template.headshape, 'color')
  if isfield(cfg.template.headshape, 'tri') && size(cfg.template.headshape.tri,1)==size(cfg.template.headshape.color,1)
    cfg.template.headshapestyle = ft_setopt(cfg.template.headshapestyle, 'facecolor', cfg.template.headshape.color);
  elseif size(cfg.template.headshape.pos,1)==size(cfg.template.headshape.color,1)
    cfg.template.headshapestyle = ft_setopt(cfg.template.headshapestyle, 'vertexcolor', cfg.template.headshape.color);
  end
end

if isempty(cfg.showalpha)
  % only show the global alpha level if not specified for any of the individual objects
  cfg.showalpha = true;
  for fn1={'individual', 'template'}
    fn1 = fn1{1};
    for fn2={'elecstyle', 'gradstyle', 'headshapestyle', 'headmodelstyle', 'mristyle', 'meshstyle'}
      fn2 = fn2{1};
      style = cfg.(fn1).(fn2);
      cfg.showalpha = cfg.showalpha & ~any(endsWith(style(1:2:end), 'alpha'));
    end
  end
end

template   = struct(cfg.template);
individual = struct(cfg.individual);

% ensure that they are consistent with the latest FieldTrip version
for fn1={'individual', 'template'}
  fn1 = fn1{1};
  for fn2={'elec', 'grad', 'opto'}
    fn2 = fn2{1};
    if ~isempty(cfg.(fn1).(fn2))
      cfg.(fn1).(fn2) = ft_datatype_sens(cfg.(fn1).(fn2));
    end
  end
end

% ensure that they are consistent with the latest FieldTrip version
if ~isempty(template.headshape)
  template.headshape = fixpos(template.headshape);
end
if ~isempty(individual.headshape)
  individual.headshape = fixpos(individual.headshape);
end

% convert the coordinates of all geometrical objects into mm
fn = {'elec', 'grad', 'headshape', 'headmodel', 'mri', 'mesh'};
hasindividual = false(size(fn));
originalunit = cell(size(fn));
for i=1:length(fn)
  if ~isempty(individual.(fn{i}))
    hasindividual(i) = true;
    originalunit{i} = individual.(fn{i}).unit;
    individual.(fn{i}) = ft_convert_units(individual.(fn{i}), cfg.unit); % ensure that the units are known and all the same
  end
end
hastemplate = false(size(fn));
for i=1:length(fn)
  if ~isempty(template.(fn{i}))
    hastemplate(i) = true;
    template.(fn{i}) = ft_convert_units(template.(fn{i}), cfg.unit); % ensure that the units are known and all the same
  end
end

% determine the coordinate system of the template objects
coordsys = [];
for i=1:length(fn)
  if ~isempty(template.(fn{i}))
    if isfield(template.(fn{i}), 'coordsys')
      if isempty(coordsys)
        % remember the first coordinate system
        coordsys = template.(fn{i}).coordsys;
      end
      % ensure that all template objects have the same coordinate system as the first one
      assert(isequal(coordsys, template.(fn{i}).coordsys));
    end
  end
end

% ensure that the headshape surface is triangulated
if ~isempty(template.headshape)
  if ~isfield(template.headshape, 'tri') || isempty(template.headshape.tri)
    template.headshape.tri = projecttri(template.headshape.pos);
  end
end
if ~isempty(individual.headshape)
  if ~isfield(individual.headshape, 'tri') || isempty(individual.headshape.tri)
    individual.headshape.tri = projecttri(individual.headshape.pos);
  end
end

ft_info('Use the mouse to rotate the geometry, and click "redisplay" to update the light.');
ft_info('Close the figure when you are done.');

% open a figure
fig = figure;
set(fig, 'CloseRequestFcn',    @cb_quit);
set(fig, 'windowkeypressfcn',  @cb_keyboard);

% move the axes so that they don't overlap with the GUI elements
set(gca, 'position', [0.15 0.30 0.50 0.60]);
% add the data and the settings to the figure
setappdata(fig, 'individual',     individual);
setappdata(fig, 'template',       template);
setappdata(fig, 'transform',      eye(4));
setappdata(fig, 'init',           true);
setappdata(fig, 'cleanup',        false);
setappdata(fig, 'rotate',         cfg.rotate);
setappdata(fig, 'scale',          cfg.scale);
setappdata(fig, 'translate',      cfg.translate);
setappdata(fig, 'coordsys',       coordsys); % can be unknown
setappdata(fig, 'unit',           cfg.unit);
setappdata(fig, 'labels',         true);
setappdata(fig, 'axes',           true);
setappdata(fig, 'grid',           true);
setappdata(fig, 'camlight',       true);
setappdata(fig, 'cfg',            cfg);

if isempty(coordsys) || strcmp(coordsys, 'unknown')
  ft_notice('the template coordinate system is unknown, selecting the viewpoint is not possible');
  ft_uilayout(gcf, 'tag', 'viewpointbtn', 'Visible', 'off');
else
  [labelx, labely, labelz] = coordsys2label(coordsys, 2, 0);
  ft_notice('the template coordinate system is "%s"', coordsys);
  ft_info('the positive X-axis is pointing to %s', labelx);
  ft_info('the positive Y-axis is pointing to %s', labely);
  ft_info('the positive Z-axis is pointing to %s', labelz);
end

% add the GUI elements
axmax = 150 * ft_scalingfactor('mm', cfg.unit);
axis([-axmax axmax -axmax axmax -axmax axmax]);
cb_creategui(gcf);
cb_redraw(gcf);
cb_help(cfg);
rotate3d on

cleanup = false;
while ~cleanup
  uiwait(fig);
  cfg.m   = getappdata(fig, 'transform');
  cleanup = getappdata(fig, 'cleanup');
end

% remember the details of the transformation
cfg.m         = getappdata(fig, 'transform');
cfg.rotate    = getappdata(fig, 'rotate');
cfg.scale     = getappdata(fig, 'scale');
cfg.translate = getappdata(fig, 'translate');

if istrue(cfg.showapply)
  % do not return the individual steps if the apply button was available and potentially pressed
  % in that case they only represent the last iteration, not the whole transformation
  cfg = removefields(cfg, {'rotate', 'scale', 'translate'});
end

delete(fig);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble provenance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_creategui(h, eventdata, handles)

% define the position of each GUI element
fig = getparent(h);

% get the initial values
rotate    = getappdata(fig, 'rotate');
scale     = getappdata(fig, 'scale');
translate = getappdata(fig, 'translate');
cfg       = getappdata(fig, 'cfg');

% constants
CONTROL_WIDTH   = 0.04;
CONTROL_HEIGHT  = 0.05;
CONTROL_HOFFSET = 0.75;
CONTROL_VOFFSET = 0.80;

% rotateui
uicontrol('tag', 'rotateui', 'parent',  fig, 'units', 'normalized', 'style', 'text', 'string', 'rotate',  'callback', [])
uicontrol('tag', 'rx',       'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', rotate(1), 'callback', @cb_redraw)
uicontrol('tag', 'ry',       'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', rotate(2), 'callback', @cb_redraw)
uicontrol('tag', 'rz',       'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', rotate(3), 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'rotateui', 'width',  2*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET,                 'vpos',  CONTROL_VOFFSET+3*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'rx',       'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+3*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'ry',       'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+3*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'rz',       'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+3*CONTROL_HEIGHT);

% scaleui
uicontrol('tag', 'scaleui', 'parent',  fig, 'units', 'normalized', 'style', 'text', 'string', 'scale',  'callback', [])
uicontrol('tag', 'sx',      'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', scale(1), 'callback', @cb_redraw)
uicontrol('tag', 'sy',      'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', scale(2), 'callback', @cb_redraw)
uicontrol('tag', 'sz',      'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', scale(3), 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'scaleui', 'width',  2*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET,                 'vpos',  CONTROL_VOFFSET+2*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sx',      'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+2*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sy',      'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+2*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sz',      'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+2*CONTROL_HEIGHT);

% translateui
uicontrol('tag', 'translateui', 'parent',  fig, 'units', 'normalized', 'style', 'text', 'string', 'translate',  'callback', [])
uicontrol('tag', 'tx',          'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', translate(1), 'callback', @cb_redraw)
uicontrol('tag', 'ty',          'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', translate(2), 'callback', @cb_redraw)
uicontrol('tag', 'tz',          'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', translate(3), 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'translateui', 'width',  2*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET,                 'vpos',  CONTROL_VOFFSET+1*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'tx',          'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+1*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'ty',          'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+1*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'tz',          'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+1*CONTROL_HEIGHT);

% alpha ui
uicontrol('tag', 'alphaui', 'parent',  fig, 'units', 'normalized', 'style', 'text', 'string', 'alpha', 'value', [], 'callback', [], 'Visible', istrue(cfg.showalpha));
uicontrol('tag', 'alpha',   'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '0.6',   'value', [], 'callback', @cb_redraw, 'Visible', istrue(cfg.showalpha));
ft_uilayout(fig, 'tag', 'alphaui', 'width',  3*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET,                 'vpos',  CONTROL_VOFFSET-1*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'alpha',   'width',  3*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET-1*CONTROL_HEIGHT);

% control GUI elements on the right side
uicontrol('tag', 'viewpointbtn',  'parent',  fig, 'units', 'normalized', 'style', 'popup',      'string', 'default|top|bottom|left|right|front|back', 'value', 1,                   'callback', @cb_viewpoint);
uicontrol('tag', 'camlightbtn',   'parent',  fig, 'units', 'normalized', 'style', 'popup',      'string', 'none|camlight|uniform',            'value', 1,                           'callback', @cb_camlight, 'Visible', istrue(cfg.showlight));
uicontrol('tag', 'materialbtn',   'parent',  fig, 'units', 'normalized', 'style', 'popup',      'string', 'shiny|dull|metal',                 'value', 1,                           'callback', @cb_material, 'Visible', istrue(cfg.showmaterial));
uicontrol('tag', 'axes1btn',      'parent',  fig, 'units', 'normalized', 'style', 'checkbox',   'string', '3d axes',                          'value', istrue(cfg.template.axes),   'callback', @cb_axes1);
uicontrol('tag', 'axes2btn',      'parent',  fig, 'units', 'normalized', 'style', 'checkbox',   'string', 'figure axes',                      'value', getappdata(fig, 'axes'),     'callback', @cb_axes2);
uicontrol('tag', 'labelsbtn',     'parent',  fig, 'units', 'normalized', 'style', 'checkbox',   'string', 'figure axes label',                'value', getappdata(fig, 'labels'),   'callback', @cb_labels);
uicontrol('tag', 'gridbtn',       'parent',  fig, 'units', 'normalized', 'style', 'checkbox',   'string', 'figure axes grid',                 'value', getappdata(fig, 'grid'),     'callback', @cb_grid);
uicontrol('tag', 'redisplaybtn',  'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'redisplay',                        'value', [],                          'callback', @cb_redisplay);
uicontrol('tag', 'applybtn',      'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'apply',                            'value', [],                          'callback', @cb_apply,    'Visible', istrue(cfg.showapply));
uicontrol('tag', 'quitbtn',       'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'quit',                             'value', 1,                           'callback', @cb_quit);

% remove the invisible GUI elements
delete(findall(fig, 'type', 'uicontrol', 'visible', false))

% organize the remaining ones
ft_uilayout(fig, 'tag', '.*btn', 'units', 'normalized', 'width', 0.25, 'height', 0.05, 'hpos', 0.75, 'vpos', 'auto', 'backgroundcolor', [0.8 0.8 0.8]);
ft_uilayout(fig, 'tag', '.*btn', 'vshift', -5*CONTROL_HEIGHT);
ft_uilayout(fig, 'style', 'checkbox', 'backgroundcolor', get(fig, 'color'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata, handles)

fig            = getparent(h);
individual     = getappdata(fig, 'individual');
template       = getappdata(fig, 'template');
transform      = getappdata(fig, 'transform');
coordsys       = getappdata(fig, 'coordsys');
unit           = getappdata(fig, 'unit');
cfg            = getappdata(fig, 'cfg');

% disable all GUI elements
set(findall(fig, 'type', 'UIControl'), 'Enable', 'off')

% get the transformation details
rx = str2double(get(findobj(fig, 'tag', 'rx'), 'string'));
ry = str2double(get(findobj(fig, 'tag', 'ry'), 'string'));
rz = str2double(get(findobj(fig, 'tag', 'rz'), 'string'));
sx = str2double(get(findobj(fig, 'tag', 'sx'), 'string'));
sy = str2double(get(findobj(fig, 'tag', 'sy'), 'string'));
sz = str2double(get(findobj(fig, 'tag', 'sz'), 'string'));
tx = str2double(get(findobj(fig, 'tag', 'tx'), 'string'));
ty = str2double(get(findobj(fig, 'tag', 'ty'), 'string'));
tz = str2double(get(findobj(fig, 'tag', 'tz'), 'string'));

% create the transformation matrix
R = rotate   ([rx ry rz]);
S = scale    ([sx sy sz]);
T = translate([tx ty tz]);
% combine the transformation from the GUI with the one that has been previously applied
H = combine_transform(R, S, T, cfg.transformorder);
transform = H * transform;

axis vis3d
cla
hold on

if getappdata(fig, 'init')
  % start from a fresh view
  [az, el] = view(3);
  setappdata(fig, 'init', false);
else
  % remember the current view
  [az, el] = view();
  % FIXME it would be possible to move the plotting from all template objects here
  % FIXME so that they do not have to be plotted again upon each update
  % FIXME however, that requires that the "cla" above is implemented more specifically
end

% the "individual" struct is a local copy, so it is safe to change it here
if ~isempty(individual.headmodel)
  individual.headmodel = ft_transform_geometry(transform, individual.headmodel);
end
if ~isempty(individual.elec)
  individual.elec = ft_transform_geometry(transform, individual.elec);
end
if ~isempty(individual.grad)
  individual.grad = ft_transform_geometry(transform, individual.grad);
end
if ~isempty(individual.headshape)
  individual.headshape = ft_transform_geometry(transform, individual.headshape);
end
if ~isempty(individual.mri)
  individual.mri = ft_transform_geometry(transform, individual.mri);
end
if ~isempty(individual.mesh)
  individual.mesh = ft_transform_geometry(transform, individual.mesh);
end

% plot all the template and individual objects
if ~isempty(template.mri)
  ft_plot_ortho(template.mri.anatomy, 'transform', template.mri.transform, 'unit', template.mri.unit, 'style', 'intersect', 'intersectmesh', individual.headshape, template.mristyle{:});
end

if ~isempty(individual.mri)
  ft_plot_ortho(individual.mri.anatomy, 'transform', individual.mri.transform, 'unit', individual.mri.unit, 'style', 'intersect', 'intersectmesh', template.headshape, individual.mristyle{:});
end

if ~isempty(template.elec)
  if isfield(template.elec, 'line')
    tmpbnd = [];
    tmpbnd.pos = template.elec.chanpos;
    tmpbnd.tri = template.elec.line;
    ft_plot_mesh(tmpbnd, 'vertexcolor', 'b', 'facecolor', 'none', 'edgecolor', 'b', 'vertexsize', 10)
  else
    ft_plot_sens(template.elec, template.elecstyle{:});
  end
end

if ~isempty(individual.elec)
  if isfield(individual.elec, 'line')
    tmpbnd = [];
    tmpbnd.pos = individual.elec.chanpos;
    tmpbnd.tri = individual.elec.line;
    ft_plot_mesh(tmpbnd, 'vertexcolor', 'r', 'facecolor', 'none', 'edgecolor', 'r', 'vertexsize', 10)
  else
    ft_plot_sens(individual.elec, individual.elecstyle{:});
  end
end

if ~isempty(template.grad)
  ft_plot_sens(template.grad, template.gradstyle{:});
end

if ~isempty(individual.grad)
  ft_plot_sens(individual.grad, individual.gradstyle{:});
end

if isstruct(template.headmodel)
  ft_plot_headmodel(template.headmodel, template.headmodelstyle{:});
end

if isstruct(individual.headmodel)
  ft_plot_headmodel(individual.headmodel, individual.headmodelstyle{:});
end

if isstruct(template.headshape) && isfield(template.headshape, 'pos') && ~isempty(template.headshape.pos)
  ft_plot_headshape(template.headshape, template.headshapestyle{:});
end

if isstruct(individual.headshape) && isfield(individual.headshape, 'pos') && ~isempty(individual.headshape.pos)
  ft_plot_headshape(individual.headshape, individual.headshapestyle{:})
end

if isstruct(template.mesh) && isfield(template.mesh, 'pos') && ~isempty(template.mesh.pos)
  % there can be multiple meshes as a struct-array
  ft_plot_mesh(template.mesh, template.meshstyle{:});
end

if isstruct(individual.mesh) && isfield(individual.mesh, 'pos') && ~isempty(individual.mesh.pos)
  % there can be multiple meshes as a struct-array
  ft_plot_mesh(individual.mesh, individual.meshstyle{:});
end

if istrue(cfg.showalpha)
  % only apply the global alpha level when it is not explicitly set for the individual objects
  alpha(str2double(get(findobj(fig, 'tag', 'alpha'), 'string')));
end

% update the figure based on the GUI elements
cb_axes1(h, []);
cb_axes2(h, []);
cb_labels(h, []);
cb_grid(h, []);
cb_camlight(h, []);
cb_material(h, []);

% restore the current view
view(az, el);

% re-enable all GUI elements
set(findall(fig, 'type', 'UIControl'), 'Enable', 'on')

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
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_material(h, eventdata)
fig = getparent(h);
switch get(findall(fig, 'tag', 'materialbtn'), 'value')
  case 1
    material shiny
  case 2
    material dull
  case 3
    material metal
end
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_axes1(h, eventdata)
fig = getparent(h);
if get(findall(fig, 'tag', 'axes1btn'), 'value')
  % FT_PLOT_AXES calls FT_PLOT_MESH, which does "axis off"
  ax = findall(fig, 'type', 'axes');
  prevaxis = ax.Visible;
  ft_plot_axes([], 'unit', getappdata(fig, 'unit'), 'coordsys', getappdata(fig, 'coordsys'), 'tag', 'axes1');
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
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_labels(h, eventdata)
fig  = getparent(h);
unit = getappdata(fig, 'unit');
if get(findall(fig, 'tag', 'labelsbtn'), 'value')
  xlabel(sprintf('x (%s)', unit))
  ylabel(sprintf('y (%s)', unit))
  zlabel(sprintf('z (%s)', unit))
else
  xlabel('')
  ylabel('')
  zlabel('')
end
uiresume;

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
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_apply(h, eventdata, handles)

fig       = getparent(h);
transform = getappdata(fig, 'transform');
cfg       = getappdata(fig, 'cfg');

% get the transformation details
rx = str2double(get(findobj(fig, 'tag', 'rx'), 'string'));
ry = str2double(get(findobj(fig, 'tag', 'ry'), 'string'));
rz = str2double(get(findobj(fig, 'tag', 'rz'), 'string'));
sx = str2double(get(findobj(fig, 'tag', 'sx'), 'string'));
sy = str2double(get(findobj(fig, 'tag', 'sy'), 'string'));
sz = str2double(get(findobj(fig, 'tag', 'sz'), 'string'));
tx = str2double(get(findobj(fig, 'tag', 'tx'), 'string'));
ty = str2double(get(findobj(fig, 'tag', 'ty'), 'string'));
tz = str2double(get(findobj(fig, 'tag', 'tz'), 'string'));

% create the transformation matrix
R = rotate   ([rx ry rz]);
S = scale    ([sx sy sz]);
T = translate([tx ty tz]);
% combine the transformation from the GUI with the one that has been previously applied
H = combine_transform(R, S, T, cfg.transformorder);
transform = H * transform;

set(findobj(fig, 'tag', 'rx'), 'string',  0);
set(findobj(fig, 'tag', 'ry'), 'string',  0);
set(findobj(fig, 'tag', 'rz'), 'string',  0);
set(findobj(fig, 'tag', 'sx'), 'string',  1);
set(findobj(fig, 'tag', 'sy'), 'string',  1);
set(findobj(fig, 'tag', 'sz'), 'string',  1);
set(findobj(fig, 'tag', 'tx'), 'string',  0);
set(findobj(fig, 'tag', 'ty'), 'string',  0);
set(findobj(fig, 'tag', 'tz'), 'string',  0);

setappdata(fig, 'transform',  transform);
setappdata(fig, 'rotate',     [rx ry rz]);
setappdata(fig, 'scale',      [sx sy sz]);
setappdata(fig, 'translate',  [tx ty tz]);

if ~getappdata(fig, 'cleanup')
  cb_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_help(h, eventdata)

fprintf('\n');
fprintf('==================================================================================\n');
fprintf('Press "h" to show this help.\n');
fprintf('Press "q" or close the window when you are done.\n');
fprintf('Press "v" to update the light position.\n');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redisplay(h, eventdata)
fig = getparent(h);
cb_viewpoint(findall(fig, 'tag', 'viewpointbtn'));
cb_redraw(h);

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
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quit(h, eventdata)
fig = getparent(h);
setappdata(fig, 'cleanup', true); % ensure that it will apply the current transformation
cb_apply(h);
uiresume;

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
function style = updatestyle(style)
if ischar(style)
  switch style
    case 'vertex'
      style = {'vertexcolor', 'k', 'edgecolor', 'none', 'facecolor', 'none'};
    case 'edge'
      style = {'vertexcolor', 'none', 'edgecolor', 'k', 'facecolor', 'none'};
    case 'surface'
      style = {'vertexcolor', 'none', 'edgecolor', 'none', 'facecolor', 'skin'};
    case 'both'
      style = {'vertexcolor', 'none', 'edgecolor', 'k', 'facecolor', 'skin'};
    otherwise
      ft_error('unsupported style "%s"', style);
  end % switch
end % if
