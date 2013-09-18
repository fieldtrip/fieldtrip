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
%  cfg.individual.vol
%  cfg.individual.elec
%  cfg.individual.grad
%  cfg.individual.headshape
%  cfg.individual.headshapestyle = 'vertex'  (default), 'surface' or 'both'
%  cfg.individual.volstyle       = 'edge'    (default), 'surface' or 'both'
%
% The configuration structure should also contain the geometrical
% objects of a template that serves as target
%  cfg.template.vol
%  cfg.template.elec
%  cfg.template.grad
%  cfg.template.headshape
%  cfg.template.headshapestyle   = 'surface' (default), 'vertex' or 'both'
%  cfg.individual.volstyle       = 'surface' (default), 'edge'   or 'both'
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', {'individual', 'template'});

if ~isfield(cfg.individual, 'vol'),              cfg.individual.vol = [];                   end
if ~isfield(cfg.individual, 'elec'),             cfg.individual.elec = [];                  end
if ~isfield(cfg.individual, 'grad'),             cfg.individual.grad = [];                  end
if ~isfield(cfg.individual, 'headshape'),        cfg.individual.headshape = [];             end
if ~isfield(cfg.individual, 'headshapestyle'),   cfg.individual.headshapestyle = 'vertex';  end
if ~isfield(cfg.individual, 'volstyle'),         cfg.individual.volstyle = 'edge';          end

if ~isfield(cfg.template, 'vol'),              cfg.template.vol = [];                     end
if ~isfield(cfg.template, 'elec'),             cfg.template.elec = [];                    end
if ~isfield(cfg.template, 'grad'),             cfg.template.grad = [];                    end
if ~isfield(cfg.template, 'headshape'),        cfg.template.headshape = [];               end
if ~isfield(cfg.template, 'headshapestyle'),   cfg.template.headshapestyle = 'surface';   end
if ~isfield(cfg.template, 'volstyle'),         cfg.template.volstyle = 'surface';         end

template   = cfg.template;
individual = cfg.individual;

if ~isempty(template.headshape)
  if ~isfield(template.headshape, 'tri') || isempty(template.headshape.tri)
    template.headshape.tri = projecttri(template.headshape.pnt);
  end
end

if ~isempty(individual.headshape) && isfield(individual.headshape, 'pnt') && ...
    ~isempty(individual.headshape.pnt)
  if ~isfield(individual.headshape, 'tri') || isempty(individual.headshape.tri)
    individual.headshape.tri = projecttri(individual.headshape.pnt);
  end
end

% open a figure
fig = figure;
% add the data to the figure
set(fig, 'CloseRequestFcn', @cb_close);
setappdata(fig, 'individual', individual);
setappdata(fig, 'template',   template);
setappdata(fig, 'transform',  eye(4));

% add the GUI elements
cb_creategui(gca);
cb_redraw(gca);
rotate3d on
waitfor(fig);

% get the data from the figure that was left behind as global variable
% FIXME pass this as appdata
global norm
tmp = norm;
clear global norm
norm = tmp;
clear tmp

% remember the transform
cfg.m = norm.m;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some simple SUBFUNCTIONs that facilitate 3D plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = my_plot3(xyz, varargin)
h = plot3(xyz(:,1), xyz(:,2), xyz(:,3), varargin{:});
function h = my_text3(xyz, varargin)
h = text(xyz(:,1), xyz(:,2), xyz(:,3), varargin{:});
function my_line3(xyzB, xyzE, varargin)
for i=1:size(xyzB,1)
  line([xyzB(i,1) xyzE(i,1)], [xyzB(i,2) xyzE(i,2)], [xyzB(i,3) xyzE(i,3)], varargin{:})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_creategui(hObject, eventdata, handles);
% define the position of each GUI element
fig = get(hObject, 'parent');
% constants
CONTROL_WIDTH  = 0.05;
CONTROL_HEIGHT = 0.06; 
CONTROL_HOFFSET = 0.7;
CONTROL_VOFFSET = 0.5;
% rotateui
uicontrol('tag', 'rotateui', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'rotate', 'callback', [])
uicontrol('tag', 'rx', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'ry', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'rz', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'rotateui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET);
ft_uilayout(fig, 'tag', 'rx',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);
ft_uilayout(fig, 'tag', 'ry',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);
ft_uilayout(fig, 'tag', 'rz',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);

% scaleui
uicontrol('tag', 'scaleui', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'scale', 'callback', [])
uicontrol('tag', 'sx', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
uicontrol('tag', 'sy', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
uicontrol('tag', 'sz', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'scaleui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sx',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sy',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'sz',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);

% translateui
uicontrol('tag', 'translateui', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'translate', 'callback', [])
uicontrol('tag', 'tx', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'ty', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'tz', 'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
ft_uilayout(fig, 'tag', 'translateui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'tx',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'ty',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(fig, 'tag', 'tz',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);

% control buttons
uicontrol('tag', 'redisplaybtn', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'redisplay', 'value', [], 'callback', @cb_redraw);
uicontrol('tag', 'applybtn', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'apply',     'value', [], 'callback', @cb_apply);
uicontrol('tag', 'toggle labels', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'toggle label', 'value', 0, 'callback', @cb_redraw);
uicontrol('tag', 'toggle axes', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'toggle axes', 'value', 0, 'callback', @cb_redraw);
ft_uilayout(fig, 'tag', 'redisplaybtn',  'BackgroundColor', [0.8 0.8 0.8], 'width', 6*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-3*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'applybtn',      'BackgroundColor', [0.8 0.8 0.8], 'width', 6*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-4*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'toggle labels', 'BackgroundColor', [0.8 0.8 0.8], 'width', 6*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-5*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'toggle axes',   'BackgroundColor', [0.8 0.8 0.8], 'width', 6*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-6*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);

% alpha ui (somehow not implemented, facealpha is fixed at 0.7
uicontrol('tag', 'alphaui', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'alpha', 'value', [], 'callback', []);
uicontrol('tag', 'alpha',   'parent', fig, 'units', 'normalized', 'style', 'edit', 'string', '0.7',   'value', [], 'callback', @cb_redraw);
ft_uilayout(fig, 'tag', 'alphaui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 3*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-7*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
ft_uilayout(fig, 'tag', 'alpha',   'BackgroundColor', [0.8 0.8 0.8], 'width', 3*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-7*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(hObject, eventdata, handles)
fig = get(hObject, 'parent');
individual = getappdata(fig, 'individual');
template = getappdata(fig, 'template');
% get the transformation details
rx = str2num(get(findobj(fig, 'tag', 'rx'), 'string'));
ry = str2num(get(findobj(fig, 'tag', 'ry'), 'string'));
rz = str2num(get(findobj(fig, 'tag', 'rz'), 'string'));
tx = str2num(get(findobj(fig, 'tag', 'tx'), 'string'));
ty = str2num(get(findobj(fig, 'tag', 'ty'), 'string'));
tz = str2num(get(findobj(fig, 'tag', 'tz'), 'string'));
sx = str2num(get(findobj(fig, 'tag', 'sx'), 'string'));
sy = str2num(get(findobj(fig, 'tag', 'sy'), 'string'));
sz = str2num(get(findobj(fig, 'tag', 'sz'), 'string'));
R = rotate   ([rx ry rz]);
T = translate([tx ty tz]);
S = scale    ([sx sy sz]);
H = S * T * R;
axis vis3d; cla
xlabel('x')
ylabel('y')
zlabel('z')

hold on

% the "individual" struct is a local copy, so it is safe to change it here
if ~isempty(individual.vol)
  individual.vol = ft_transform_vol(H, individual.vol);
end
if ~isempty(individual.elec)
  individual.elec = ft_transform_sens(H, individual.elec);
end
if ~isempty(individual.grad)
  individual.grad = ft_transform_sens(H, individual.grad);
end
if ~isempty(individual.headshape)
  individual.headshape = ft_transform_headshape(H, individual.headshape);
end

if ~isempty(template.elec)
  if isfield(template.elec, 'line')
    tmpbnd = [];
    tmpbnd.pnt = template.elec.chanpos;
    tmpbnd.tri = template.elec.line;
    ft_plot_mesh(tmpbnd,'vertexcolor','b','facecolor','none','edgecolor','b','vertexsize',10)
  else
    ft_plot_mesh(template.elec.chanpos,'vertexcolor','b','vertexsize',10);
  end
end

if ~isempty(individual.elec)
  if isfield(individual.elec, 'line')
    tmpbnd = [];
    tmpbnd.pnt = individual.elec.chanpos;
    tmpbnd.tri = individual.elec.line;
    ft_plot_mesh(tmpbnd,'vertexcolor','r','facecolor','none','edgecolor','r','vertexsize',10)
  else
    ft_plot_mesh(individual.elec.chanpos,'vertexcolor','r','vertexsize',10);
  end
end

if ~isempty(template.grad)
  ft_plot_mesh(template.grad.chanpos,'vertexcolor','b','vertexsize',10);
  % FIXME also plot lines?
end

if ~isempty(individual.grad)
  ft_plot_mesh(individual.grad.chanpos,'vertexcolor','r','vertexsize',10);
  % FIXME also plot lines?
end

if ~isempty(template.vol)
  % FIXME this only works for boundary element models
  for i = 1:numel(template.vol.bnd)
    if strcmp(template.volstyle, 'edge') || ...
        strcmp(template.volstyle, 'both')
      ft_plot_mesh(template.vol.bnd(i),'facecolor','none','vertexcolor','b')
    end
    if strcmp(template.volstyle, 'surface') || ...
        strcmp(template.volstyle, 'both')
      ft_plot_mesh(template.vol.bnd(i),'facecolor','b','edgecolor','none')
      lighting gouraud
      material shiny
      camlight
    end
  end
end

if ~isempty(individual.vol)
  % FIXME this only works for boundary element models
  for i = 1:numel(individual.vol.bnd)
    if strcmp(individual.volstyle, 'edge') || ...
        strcmp(individual.volstyle, 'both')
      ft_plot_mesh(individual.vol.bnd(i),'facecolor','none','vertexcolor','r')
    end
    if strcmp(individual.volstyle, 'surface') || ...
        strcmp(individual.volstyle, 'both')
      ft_plot_mesh(individual.vol.bnd(i),'facecolor','r','edgecolor','none')
      lighting gouraud
      material shiny
      camlight
    end
  end
end

if ~isempty(template.headshape)
  if isfield(template.headshape, 'pnt') && ~isempty(template.headshape.pnt)
    if strcmp(template.headshapestyle, 'surface') || ...
        strcmp(template.headshapestyle, 'both')
      ft_plot_mesh(template.headshape,'facecolor','b','edgecolor','none')
      lighting gouraud
      material shiny
      camlight
      alpha(str2num(get(findobj(fig, 'tag', 'alpha'), 'string')));
    end
    
    if strcmp(template.headshapestyle, 'vertex') || ...
        strcmp(template.headshapestyle, 'both')
      ft_plot_mesh(template.headshape.pnt,'vertexcolor','b')
    end
  end
  if isfield(template.headshape, 'fid') && ~isempty(template.headshape.fid.pnt)
    ft_plot_mesh(template.headshape.fid.pnt,'vertexcolor','b','vertexsize',20)
  end
end

if ~isempty(individual.headshape)
  if isfield(individual.headshape, 'pnt') && ~isempty(individual.headshape.pnt)
    if strcmp(individual.headshapestyle, 'surface') || ...
        strcmp(individual.headshapestyle, 'both')
      ft_plot_mesh(individual.headshape,'facecolor','r','edgecolor','none')
      lighting gouraud
      material shiny
      camlight
      alpha(str2num(get(findobj(fig, 'tag', 'alpha'), 'string')));
    end
    
    if strcmp(individual.headshapestyle, 'vertex') || ...
        strcmp(individual.headshapestyle, 'both')
      ft_plot_mesh(individual.headshape.pnt,'vertexcolor','r')
    end
  end
  if isfield(individual.headshape, 'fid') && ~isempty(individual.headshape.fid.pnt)
    ft_plot_mesh(individual.headshape.fid.pnt,'vertexcolor','r','vertexsize',20)
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
function cb_apply(hObject, eventdata, handles);
fig = get(hObject, 'parent');
transform = getappdata(fig, 'transform');
% get the transformation details
rx = str2num(get(findobj(fig, 'tag', 'rx'), 'string'));
ry = str2num(get(findobj(fig, 'tag', 'ry'), 'string'));
rz = str2num(get(findobj(fig, 'tag', 'rz'), 'string'));
tx = str2num(get(findobj(fig, 'tag', 'tx'), 'string'));
ty = str2num(get(findobj(fig, 'tag', 'ty'), 'string'));
tz = str2num(get(findobj(fig, 'tag', 'tz'), 'string'));
sx = str2num(get(findobj(fig, 'tag', 'sx'), 'string'));
sy = str2num(get(findobj(fig, 'tag', 'sy'), 'string'));
sz = str2num(get(findobj(fig, 'tag', 'sz'), 'string'));
R = rotate   ([rx ry rz]);
T = translate([tx ty tz]);
S = scale    ([sx sy sz]);
H = S * T * R;
transform = H * transform;
set(findobj(fig, 'tag', 'rx'), 'string', 0);
set(findobj(fig, 'tag', 'ry'), 'string', 0);
set(findobj(fig, 'tag', 'rz'), 'string', 0);
set(findobj(fig, 'tag', 'tx'), 'string', 0);
set(findobj(fig, 'tag', 'ty'), 'string', 0);
set(findobj(fig, 'tag', 'tz'), 'string', 0);
set(findobj(fig, 'tag', 'sx'), 'string', 1);
set(findobj(fig, 'tag', 'sy'), 'string', 1);
set(findobj(fig, 'tag', 'sz'), 'string', 1);
setappdata(fig, 'transform', transform);
cb_redraw(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_close(hObject, eventdata, handles);
% make the current transformation permanent and subsequently allow deleting the figure
cb_apply(gca);
% get the updated electrode from the figure
fig    = hObject;
% hmmm, this is ugly
global norm
norm.m = getappdata(fig, 'transform');
set(fig, 'CloseRequestFcn', @delete);
delete(fig);

