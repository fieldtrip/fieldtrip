function cfg = ft_interactiverealign(cfg)

% FT_INTERACTIVEREALIGN interactively rotates, scales and translates
% electrode positions to template electrode positions or towards
% the head surface.
%
% Use as
%   [cfg] = ft_interactiverealign(cfg)
%
% The configuration structure should contain the individuals geometrical
% objects that have to be realigned as
%  cfg.individual.elec
%  cfg.individual.grad
%  cfg.individual.headmodel
%  cfg.individual.headmodelstyle = 'edge'    (default), 'surface' or 'both'
%  cfg.individual.headshape
%  cfg.individual.headshapestyle = 'vertex'  (default), 'surface' or 'both'
%
% The configuration structure should also contain the geometrical
% objects of a template that serves as target
%  cfg.template.elec
%  cfg.template.grad
%  cfg.template.headmodel
%  cfg.template.headmodelstyle   = 'surface' (default), 'edge'   or 'both'
%  cfg.template.headshape
%  cfg.template.headshapestyle   = 'surface' (default), 'vertex' or 'both'
%
% See also FT_VOLUMEREALIGN, FT_ELECTRODEREALIGN, FT_READ_SENS, FT_READ_VOL, FT_READ_HEADSHAPE

% Copyright (C) 2008, Vladimir Litvak
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', {'individual', 'template'});
cfg = ft_checkconfig(cfg, 'required', {'individual', 'template'});
cfg = ft_checkconfig(cfg.individual, 'renamed', {'vol', 'headmodel'});
cfg = ft_checkconfig(cfg.individual, 'renamed', {'volstyle', 'headmodelstyle'});
cfg = ft_checkconfig(cfg.template, 'renamed', {'vol', 'headmodel'});
cfg = ft_checkconfig(cfg.template, 'renamed', {'volstyle', 'headmodelstyle'});

cfg.individual.elec           = ft_getopt(cfg.individual, 'elec', []);
cfg.individual.grad           = ft_getopt(cfg.individual, 'grad', []);
cfg.individual.headshape      = ft_getopt(cfg.individual, 'headshape', []);
cfg.individual.headshapestyle = ft_getopt(cfg.individual, 'headshapestyle', 'vertex');
cfg.individual.headmodel      = ft_getopt(cfg.individual, 'headmodel', []);
cfg.individual.headmodelstyle = ft_getopt(cfg.individual, 'headmodelstyle', 'edge');
cfg.individual.mri            = ft_getopt(cfg.individual, 'mri', []);
cfg.individual.mristyle       = ft_getopt(cfg.individual, 'mristyle', 'intersect');

cfg.template.elec             = ft_getopt(cfg.template, 'elec', []);
cfg.template.grad             = ft_getopt(cfg.template, 'grad', []);
cfg.template.headshape        = ft_getopt(cfg.template, 'headshape', []);
cfg.template.headshapestyle   = ft_getopt(cfg.template, 'headshapestyle', 'vertex');
cfg.template.headmodel        = ft_getopt(cfg.template, 'headmodel', []);
cfg.template.headmodelstyle   = ft_getopt(cfg.template, 'headmodelstyle', 'edge');
cfg.template.mri              = ft_getopt(cfg.template, 'mri', []);
cfg.template.mristyle       = ft_getopt(cfg.template, 'mristyle', 'intersect');

template   = struct(cfg.template);
individual = struct(cfg.individual);

% convert the coordinates of all geometrical objects into mm
fn = {'elec', 'grad', 'headshape', 'headmodel', 'mri'};
for i=1:length(fn)
  if ~isempty(cfg.individual.(fn{i}))
    cfg.individual.(fn{i}) = ft_convert_units(cfg.individual.(fn{i}), 'mm');
  end
end
for i=1:length(fn)
  if ~isempty(cfg.template.(fn{i}))
    cfg.template.(fn{i}) = ft_convert_units(cfg.template.(fn{i}), 'mm');
  end
end

if ~isempty(template.headshape)
  if ~isfield(template.headshape, 'tri') || isempty(template.headshape.tri)
    template.headshape.tri = projecttri(template.headshape.pnt);
  end
end

if ~isempty(individual.headshape) && isfield(individual.headshape, 'pnt') && ~isempty(individual.headshape.pnt)
  if ~isfield(individual.headshape, 'tri') || isempty(individual.headshape.tri)
    individual.headshape.tri = projecttri(individual.headshape.pnt);
  end
end

% open a figure
fig = figure;
set(gca, 'position', [0.05 0.15 0.75 0.75]);

% add the data to the figure
set(fig, 'CloseRequestFcn', @cb_quit);
setappdata(fig, 'individual',  individual);
setappdata(fig, 'template',    template);
setappdata(fig, 'transform',   eye(4));
setappdata(fig, 'cleanup',     false);

% add the GUI elements
cb_creategui(gca);
cb_redraw(gca);
rotate3d on
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
uicontrol('tag', 'viewbtn',   'parent',  fig, 'units', 'normalized', 'style', 'popup', 'string', 'top|bottom|left|right|front|back', 'value',  1, 'callback', @cb_view);
uicontrol('tag', 'redisplaybtn', 'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'redisplay', 'value', [], 'callback', @cb_redraw);
uicontrol('tag', 'applybtn',  'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'apply',    'value', [], 'callback', @cb_apply);
uicontrol('tag', 'toggle labels', 'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'toggle label', 'value',  0,  'callback', @cb_redraw);
uicontrol('tag', 'toggle axes', 'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'toggle axes', 'value',  0,  'callback', @cb_redraw);
uicontrol('tag', 'quitbtn',   'parent',  fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'quit',     'value',  1,  'callback', @cb_quit);
ft_uilayout(fig, 'tag', 'viewbtn',   'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-2*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'redisplaybtn', 'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-3*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'applybtn',  'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-4*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'toggle labels', 'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-5*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'toggle axes', 'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-6*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'quitbtn',   'BackgroundColor', [0.8 0.8 0.8], 'width',  6*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-7*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);

% alpha ui (somehow not implemented, facealpha is fixed at 0.7
uicontrol('tag', 'alphaui', 'parent',  fig, 'units', 'normalized', 'style', 'text', 'string', 'alpha', 'value', [], 'callback', []);
uicontrol('tag', 'alpha', 'parent',  fig, 'units', 'normalized', 'style', 'edit', 'string', '0.6', 'value', [], 'callback', @cb_redraw);
ft_uilayout(fig, 'tag', 'alphaui', 'BackgroundColor', [0.8 0.8 0.8], 'width',  3*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-8*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'alpha', 'BackgroundColor', [0.8 0.8 0.8], 'width',  3*CONTROL_WIDTH, 'height',  CONTROL_HEIGHT, 'vpos',  CONTROL_VOFFSET-8*CONTROL_HEIGHT, 'hpos',  CONTROL_HOFFSET+3*CONTROL_WIDTH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata, handles)

fig        = getparent(h);
individual = getappdata(fig, 'individual');
template   = getappdata(fig, 'template');
transform  = getappdata(fig, 'transform');

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
xlabel('x')
ylabel('y')
zlabel('z')

hold on

% the "individual" struct is a local copy, so it is safe to change it here
if ~isempty(individual.headmodel)
  individual.headmodel = ft_transform_vol(transform, individual.headmodel);
end
if ~isempty(individual.elec)
  individual.elec = ft_transform_sens(transform, individual.elec);
end
if ~isempty(individual.grad)
  individual.grad = ft_transform_sens(transform, individual.grad);
end
if ~isempty(individual.headshape)
  individual.headshape = ft_transform_headshape(transform, individual.headshape);
end
if ~isempty(individual.mri)
  individual.mri.transform = transform * individual.mri.transform;
end

if ~isempty(template.mri)
  if strcmp(template.mristyle, 'intersect')
    ft_plot_ortho(template.mri.anatomy, 'transform',  template.mri.transform, 'style', 'intersect', 'intersectmesh',  individual.headshape);
  elseif strcmp(template.mristyle, 'montage')
    ft_plot_montage(template.mri.anatomy, 'transform',  template.mri.transform, 'style', 'intersect', 'intersectmesh',  individual.headshape);
  end
end

if ~isempty(individual.mri)
  if strcmp(individual.mristyle, 'intersect')
    ft_plot_ortho(individual.mri.anatomy, 'transform',  individual.mri.transform, 'style', 'intersect', 'intersectmesh',  template.headshape);
  elseif strcmp(individual.mristyle, 'montage')
    ft_plot_montage(individual.mri.anatomy, 'transform',  individual.mri.transform, 'style', 'intersect', 'intersectmesh',  template.headshape);
  end
end

if ~isempty(template.elec)
  % FIXME use ft_plot_sens
  if isfield(template.elec, 'line')
    tmpbnd = [];
    tmpbnd.pnt = template.elec.chanpos;
    tmpbnd.tri = template.elec.line;
    ft_plot_mesh(tmpbnd,'vertexcolor', 'b', 'facecolor', 'none', 'edgecolor', 'b', 'vertexsize',10)
  else
    ft_plot_mesh(template.elec.chanpos,'vertexcolor', 'b', 'vertexsize',10);
  end
end

if ~isempty(individual.elec)
  % FIXME use ft_plot_sens
  if isfield(individual.elec, 'line')
    tmpbnd = [];
    tmpbnd.pnt = individual.elec.chanpos;
    tmpbnd.tri = individual.elec.line;
    ft_plot_mesh(tmpbnd,'vertexcolor', 'r', 'facecolor', 'none', 'edgecolor', 'r', 'vertexsize',10)
  else
    ft_plot_mesh(individual.elec.chanpos,'vertexcolor', 'r', 'vertexsize',10);
  end
end

if ~isempty(template.grad)
  % FIXME use ft_plot_sens
  ft_plot_mesh(template.grad.chanpos,'vertexcolor', 'b', 'vertexsize',10);
  % FIXME also plot lines?
end

if ~isempty(individual.grad)
  % FIXME use ft_plot_sens
  ft_plot_mesh(individual.grad.chanpos,'vertexcolor', 'r', 'vertexsize',10);
  % FIXME also plot lines?
end

if ~isempty(template.headmodel)
  % FIXME this only works for boundary element models
  for i = 1:numel(template.headmodel.bnd)
    if strcmp(template.headmodelstyle, 'edge') || strcmp(template.headmodelstyle, 'both')
      ft_plot_mesh(template.headmodel.bnd(i), 'facecolor', 'none', 'vertexcolor', 'b')
    end
    if strcmp(template.headmodelstyle, 'surface') || ft_plot_mesh(template.headmodel.bnd(i), 'facecolor', 'b', 'edgecolor', 'none')
      lighting gouraud
      material shiny
      camlight
    end
  end
end

if ~isempty(individual.headmodel)
  % FIXME this only works for boundary element models
  for i = 1:numel(individual.headmodel.bnd)
    if strcmp(individual.headmodelstyle, 'edge') || strcmp(individual.headmodelstyle, 'both')
      ft_plot_mesh(individual.headmodel.bnd(i), 'facecolor', 'none', 'vertexcolor', 'r')
    end
    if strcmp(individual.headmodelstyle, 'surface') || strcmp(individual.headmodelstyle, 'both')
      ft_plot_mesh(individual.headmodel.bnd(i), 'facecolor', 'r', 'edgecolor', 'none')
      lighting gouraud
      material shiny
      camlight
    end
  end
end

if ~isempty(template.headshape)
  if isfield(template.headshape, 'pnt') && ~isempty(template.headshape.pnt)
    if strcmp(template.headshapestyle, 'surface') || strcmp(template.headshapestyle, 'both')
      ft_plot_mesh(template.headshape,'facecolor', 'b', 'edgecolor', 'none')
      lighting gouraud
      material shiny
      camlight
      alpha(str2double(get(findobj(fig, 'tag', 'alpha'), 'string')));
    end
    
    if strcmp(template.headshapestyle, 'vertex') || strcmp(template.headshapestyle, 'both')
      ft_plot_mesh(template.headshape.pnt,'vertexcolor', 'b')
    end
  end
  if isfield(template.headshape, 'fid') && ~isempty(template.headshape.fid.pnt)
    ft_plot_mesh(template.headshape.fid.pnt,'vertexcolor', 'b', 'vertexsize',20)
  end
end

if ~isempty(individual.headshape)
  if isfield(individual.headshape, 'pnt') && ~isempty(individual.headshape.pnt)
    if strcmp(individual.headshapestyle, 'surface') || strcmp(individual.headshapestyle, 'both')
      ft_plot_mesh(individual.headshape,'facecolor', 'r', 'edgecolor', 'none')
      lighting gouraud
      material shiny
      camlight
      alpha(str2double(get(findobj(fig, 'tag', 'alpha'), 'string')));
    end
    
    if strcmp(individual.headshapestyle, 'vertex') || strcmp(individual.headshapestyle, 'both')
      ft_plot_mesh(individual.headshape.pnt,'vertexcolor', 'r')
    end
  end
  if isfield(individual.headshape, 'fid') && ~isempty(individual.headshape.fid.pnt)
    ft_plot_mesh(individual.headshape.fid.pnt,'vertexcolor', 'r', 'vertexsize',20)
  end
end

if get(findobj(fig, 'tag', 'toggle axes'), 'value')
  axis on
else
  axis off
end
if get(findobj(fig, 'tag', 'toggle grid'), 'value')
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
function cb_view(h, eventdata)

val = get(h, 'value');
switch val
  case 1
    view([90 90]);
  case 2
    view([90 -90]);
  case 3
    view([-90 0]);
  case 4
    view([90 0]);
  case 5
    view([-180 0]);
  case 6
    view([0 0]);
    
  otherwise
end
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
