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
%   cfg.individual.elec           = structure
%   cfg.individual.grad           = structure
%   cfg.individual.headmodel      = structure, see FT_PREPARE_HEADMODEL
%   cfg.individual.headshape      = structure, see FT_READ_HEADSHAPE
%   cfg.individual.mri            = structure, see FT_READ_MRI
% You can specify the style with which the objects are displayed using
%   cfg.individual.headmodelstyle = 'vertex', 'edge', 'surface' or 'both' (default = 'edge')
%   cfg.individual.headshapestyle = 'vertex', 'edge', 'surface' or 'both' (default = 'vertex')
%
% The configuration structure should also contain the geometrical object of a
% template that serves as target
%   cfg.template.axes           = string, 'yes' or 'no (default = 'no')
%   cfg.template.elec           = structure
%   cfg.template.grad           = structure
%   cfg.template.headmodel      = structure, see FT_PREPARE_HEADMODEL
%   cfg.template.headshape      = structure, see FT_READ_HEADSHAPE
%   cfg.template.mri            = structure, see FT_READ_MRI
% You can specify the style with which the objects are displayed using
%   cfg.template.headmodelstyle = 'vertex', 'edge', 'surface' or 'both' (default = 'edge')
%   cfg.template.headshapestyle = 'vertex', 'edge', 'surface' or 'both' (default = 'vertex')
%
% You can specify one or multiple individual objects which will all be realigned and
% one or multiple template objects.
%
% See also FT_VOLUMEREALIGN, FT_ELECTRODEREALIGN, FT_DETERMINE_COORDSYS,
% FT_READ_SENS, FT_READ_VOL, FT_READ_HEADSHAPE

% Copyright (C) 2008, Vladimir Litvak
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
ft_preamble trackconfig

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

cfg.individual.elec           = ft_getopt(cfg.individual, 'elec', []);
cfg.individual.elecstyle      = ft_getopt(cfg.individual, 'elecstyle', {}); % key-value pairs
cfg.individual.grad           = ft_getopt(cfg.individual, 'grad', []);
cfg.individual.gradstyle      = ft_getopt(cfg.individual, 'gradstyle', {}); % key-value pairs
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

cfg.template.axes             = ft_getopt(cfg.template, 'axes', 'no');
cfg.template.elec             = ft_getopt(cfg.template, 'elec', []);
cfg.template.elecstyle        = ft_getopt(cfg.template, 'elecstyle', {}); % key-value pairs
cfg.template.grad             = ft_getopt(cfg.template, 'grad', []);
cfg.template.gradstyle        = ft_getopt(cfg.template, 'gradstyle', {}); % key-value pairs
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

% convert the string that describes the style to a cell-array
cfg.template.headshapestyle   = updatestyle(cfg.template.headshapestyle);
cfg.individual.headshapestyle = updatestyle(cfg.individual.headshapestyle);
cfg.template.headmodelstyle   = updatestyle(cfg.template.headmodelstyle);
cfg.individual.headmodelstyle = updatestyle(cfg.individual.headmodelstyle);

template   = struct(cfg.template);
individual = struct(cfg.individual);

% ensure that they are consistent with the latest FieldTrip version
if ~isempty(template.elec)
  template.elec = ft_datatype_sens(template.elec);
end
if ~isempty(individual.elec)
  individual.elec = ft_datatype_sens(individual.elec);
end
if ~isempty(template.headshape)
  template.headshape = fixpos(template.headshape);
end
if ~isempty(individual.headshape)
  individual.headshape = fixpos(individual.headshape);
end

% convert the coordinates of all geometrical objects into mm
fn = {'elec', 'grad', 'headshape', 'headmodel', 'mri'};
hasindividual = false(size(fn));
originalunit = cell(size(fn));
for i=1:length(fn)
  if ~isempty(individual.(fn{i}))
    hasindividual(i) = true;
    individual.(fn{i}) = ft_determine_units(individual.(fn{i})); % ensure that it has units
    originalunit{i} = individual.(fn{i}).unit;
    individual.(fn{i}) = ft_convert_units(individual.(fn{i}), 'mm'); % ensure that the units are all in mm
  end
end
hastemplate = false(size(fn));
for i=1:length(fn)
  if ~isempty(template.(fn{i}))
    template.(fn{i}) = ft_determine_units(template.(fn{i})); % ensure that it has units
    hastemplate(i) = true;
    template.(fn{i}) = ft_convert_units(template.(fn{i}), 'mm'); % ensure that the units are all in mm
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
set(gca, 'position', [0.05 0.15 0.75 0.75]);
axis([-150 150 -150 150 -150 150]);

% add the data to the figure
set(fig, 'CloseRequestFcn', @cb_quit);
setappdata(fig, 'individual',  individual);
setappdata(fig, 'template',    template);
setappdata(fig, 'transform',   eye(4));
setappdata(fig, 'cleanup',     false);
setappdata(fig, 'coordsys',    coordsys); % can be unknown

setappdata(fig, 'toggle_axes', 1);
setappdata(fig, 'toggle_grid', 1);

% add the GUI elements
cb_creategui(gcf);
cb_redraw(gcf);
rotate3d on

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

cleanup = false;
while ~cleanup
  uiwait(fig);
  cfg.m   = getappdata(fig, 'transform');
  cleanup = getappdata(fig, 'cleanup');
end

% remember the transform and touch it
cfg.m = getappdata(fig, 'transform');
cfg.m;

delete(fig);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_creategui(h, eventdata, handles)

% define the position of each GUI element
fig = getparent(h);

% constants
CONTROL_WIDTH  = 0.04;
CONTROL_HEIGHT = 0.05;
CONTROL_HOFFSET = 0.75;
CONTROL_VOFFSET = 0.5;

% rotateui
uicontrol('tag', 'rotateui', 'parent',  fig, 'units', 'normalized', 'style', 'text', 'string', 'rotate', 'callback', [])
uicontrol('tag', 'rx',   'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'ry',   'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'rz',   'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'rotateui', 'BackgroundColor', [0.8 0.8 0.8], 'width',  2*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET,                 'vpos',  CONTROL_VOFFSET+CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'rx',   'BackgroundColor', [0.8 0.8 0.8], 'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'ry',   'BackgroundColor', [0.8 0.8 0.8], 'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'rz',   'BackgroundColor', [0.8 0.8 0.8], 'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET+CONTROL_HEIGHT);

% scaleui
uicontrol('tag', 'scaleui', 'parent',  fig, 'units', 'normalized', 'style', 'text', 'string', 'scale', 'callback', [])
uicontrol('tag', 'sx',  'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
uicontrol('tag', 'sy',  'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
uicontrol('tag', 'sz',  'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'scaleui', 'BackgroundColor', [0.8 0.8 0.8], 'width',  2*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET,                 'vpos',  CONTROL_VOFFSET-0*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sx',  'BackgroundColor', [0.8 0.8 0.8], 'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET-0*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sy',  'BackgroundColor', [0.8 0.8 0.8], 'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET-0*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sz',  'BackgroundColor', [0.8 0.8 0.8], 'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET-0*CONTROL_HEIGHT);

% translateui
uicontrol('tag', 'translateui', 'parent',  fig, 'units', 'normalized', 'style', 'text', 'string', 'translate', 'callback', [])
uicontrol('tag', 'tx',      'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'ty',      'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'tz',      'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'translateui', 'BackgroundColor', [0.8 0.8 0.8], 'width',  2*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET,                 'vpos',  CONTROL_VOFFSET-1*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'tx',      'BackgroundColor', [0.8 0.8 0.8], 'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET-1*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'ty',      'BackgroundColor', [0.8 0.8 0.8], 'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET-1*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'tz',      'BackgroundColor', [0.8 0.8 0.8], 'width',  CONTROL_WIDTH,   'height',  CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos',  CONTROL_VOFFSET-1*CONTROL_HEIGHT);

% control buttons
uicontrol('tag', 'viewpointbtn',  'parent',  fig, 'units', 'normalized', 'style', 'popup',      'string', 'top|bottom|left|right|front|back', 'value',  1, 'callback', @cb_viewpoint);
uicontrol('tag', 'redisplaybtn',  'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'redisplay',    'value', [], 'callback', @cb_redraw);
uicontrol('tag', 'applybtn',      'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'apply',        'value', [], 'callback', @cb_apply);
uicontrol('tag', 'toggle labels', 'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'toggle label', 'value',  0,  'callback', @cb_redraw);
uicontrol('tag', 'toggle axes',   'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'toggle axes',  'value',  getappdata(fig, 'toggle_axes'),  'callback', @cb_redraw);
uicontrol('tag', 'toggle grid',   'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'toggle grid',  'value',  getappdata(fig, 'toggle_grid'),  'callback', @cb_redraw);
uicontrol('tag', 'quitbtn',       'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'quit',         'value',  1,  'callback', @cb_quit);
ft_uilayout(fig, 'tag', 'viewpointbtn',   'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-2*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'redisplaybtn',   'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-4*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'applybtn',       'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-5*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'toggle labels',  'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-6*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'toggle axes',    'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-7*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'toggle grid',    'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-8*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'quitbtn',        'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-9*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);

% alpha ui (somehow not implemented, facealpha is fixed at 0.7
uicontrol('tag', 'alphaui', 'parent',  fig, 'units', 'normalized', 'style', 'text', 'string', 'alpha', 'value', [], 'callback', []);
uicontrol('tag', 'alpha',   'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '0.6', 'value', [], 'callback', @cb_redraw);
ft_uilayout(fig, 'tag', 'alphaui',  'BackgroundColor', [0.8 0.8 0.8], 'width',  3*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-3*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'alpha',    'BackgroundColor', [0.8 0.8 0.8], 'width',  3*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-3*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+3*CONTROL_WIDTH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata, handles)

fig        = getparent(h);
individual = getappdata(fig, 'individual');
template   = getappdata(fig, 'template');
transform  = getappdata(fig, 'transform');
coordsys   = getappdata(fig, 'coordsys');

% get the transformation details
rx = str2double(get(findobj(fig, 'tag', 'rx'), 'string'));
ry = str2double(get(findobj(fig, 'tag', 'ry'), 'string'));
rz = str2double(get(findobj(fig, 'tag', 'rz'), 'string'));
tx = str2double(get(findobj(fig, 'tag', 'tx'), 'string'));
ty = str2double(get(findobj(fig, 'tag', 'ty'), 'string'));
tz = str2double(get(findobj(fig, 'tag', 'tz'), 'string'));
sx = str2double(get(findobj(fig, 'tag', 'sx'), 'string'));
sy = str2double(get(findobj(fig, 'tag', 'sy'), 'string'));
sz = str2double(get(findobj(fig, 'tag', 'sz'), 'string'));

R = rotate   ([rx ry rz]);
T = translate([tx ty tz]);
S = scale    ([sx sy sz]);
H = S * T * R;
% combine the present transform according to the GUI with the one that has been previously applied
transform = H * transform;

axis vis3d; cla
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

hold on

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

if ~isempty(template.mri)
  ft_plot_ortho(template.mri.anatomy, 'transform',  template.mri.transform, 'style', 'intersect', 'intersectmesh', individual.headshape, cfg.individual.mristyle{:});
end

if ~isempty(individual.mri)
  ft_plot_ortho(individual.mri.anatomy, 'transform',  individual.mri.transform, 'style', 'intersect', 'intersectmesh', template.headshape, cfg.template.mristyle{:});
end

if istrue(template.axes)
  ft_plot_axes([], 'unit', 'mm', 'coordsys', coordsys);
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

% FIXME this only works for boundary element models
if isstruct(template.headmodel) && isfield(template.headmodel, 'bnd')
  for i = 1:numel(template.headmodel.bnd)
    ft_plot_mesh(template.headmodel.bnd(i), template.headmodelstyle{:});
  end
end

% FIXME this only works for boundary element models
if isstruct(individual.headmodel) && isfield(individual.headmodel, 'bnd')
  for i = 1:numel(individual.headmodel.bnd)
    ft_plot_mesh(individual.headmodel.bnd(i), individual.headmodelstyle{:});
  end
end

if isstruct(template.headshape) && isfield(template.headshape, 'pos') && ~isempty(template.headshape.pos)
  ft_plot_headshape(template.headshape, template.headshapestyle{:});
end

if isstruct(individual.headshape) && isfield(individual.headshape, 'pos') && ~isempty(individual.headshape.pos)
  ft_plot_headshape(individual.headshape, individual.headshapestyle{:})
end

alpha(str2double(get(findobj(fig, 'tag', 'alpha'), 'string')));
lighting gouraud
material shiny
camlight

if strcmp(get(h, 'tag'), 'toggle axes')
  setappdata(fig, 'toggle_axes', ~getappdata(fig, 'toggle_axes'))
end

if getappdata(fig, 'toggle_axes')
  axis on
else
  axis off
end

if strcmp(get(h, 'tag'), 'toggle grid')
  setappdata(fig, 'toggle_grid', ~getappdata(fig, 'toggle_grid'))
end

if getappdata(fig, 'toggle_grid')
  grid on
else
  grid off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_apply(h, eventdata, handles)

fig       = getparent(h);
transform = getappdata(fig, 'transform');

% get the transformation details
rx = str2double(get(findobj(fig, 'tag', 'rx'), 'string'));
ry = str2double(get(findobj(fig, 'tag', 'ry'), 'string'));
rz = str2double(get(findobj(fig, 'tag', 'rz'), 'string'));
tx = str2double(get(findobj(fig, 'tag', 'tx'), 'string'));
ty = str2double(get(findobj(fig, 'tag', 'ty'), 'string'));
tz = str2double(get(findobj(fig, 'tag', 'tz'), 'string'));
sx = str2double(get(findobj(fig, 'tag', 'sx'), 'string'));
sy = str2double(get(findobj(fig, 'tag', 'sy'), 'string'));
sz = str2double(get(findobj(fig, 'tag', 'sz'), 'string'));

% create the transformation matrix;
R = rotate   ([rx ry rz]);
T = translate([tx ty tz]);
S = scale    ([sx sy sz]);
H = S * T * R;
transform = H * transform;

set(findobj(fig, 'tag', 'rx'), 'string',  0);
set(findobj(fig, 'tag', 'ry'), 'string',  0);
set(findobj(fig, 'tag', 'rz'), 'string',  0);
set(findobj(fig, 'tag', 'tx'), 'string',  0);
set(findobj(fig, 'tag', 'ty'), 'string',  0);
set(findobj(fig, 'tag', 'tz'), 'string',  0);
set(findobj(fig, 'tag', 'sx'), 'string',  1);
set(findobj(fig, 'tag', 'sy'), 'string',  1);
set(findobj(fig, 'tag', 'sz'), 'string',  1);

setappdata(fig, 'transform',  transform);

if ~getappdata(fig, 'cleanup')
  cb_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_viewpoint(h, eventdata)

fig       = getparent(h);
coordsys  = getappdata(fig, 'coordsys');

% get the index of the option that was selected
val = get(h, 'value');

viewpoint = {'top', 'bottom', 'left', 'right', 'front', 'back'};
setviewpoint(gca, coordsys, viewpoint{val});

uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quit(h, eventdata)

fig = getparent(h);
setappdata(fig, 'cleanup',  true);

% ensure to apply the current transformation
cb_apply(h);

uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

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
