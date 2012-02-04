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
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

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
ft_postamble trackconfig
ft_postamble callinfo


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
% SUBFUNCTION to layout a moderately complex graphical user interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = layoutgui(fig, geometry, position, style, string, value, tag, callback);
horipos  = geometry(1); % lower left corner of the GUI part in the figure
vertpos  = geometry(2); % lower left corner of the GUI part in the figure
width    = geometry(3); % width  of the GUI part in the figure
height   = geometry(4); % height of the GUI part in the figure
horidist = 0.05;
vertdist = 0.05;
options  = {'units', 'normalized', 'HorizontalAlignment', 'center'}; %  'VerticalAlignment', 'middle'
Nrow     = size(position,1);
h        = cell(Nrow,1);
for i=1:Nrow
  if isempty(position{i})
    continue;
  end
  position{i} = position{i} ./ sum(position{i});
  Ncol = size(position{i},2);
  ybeg = (Nrow-i  )/Nrow + vertdist/2;
  yend = (Nrow-i+1)/Nrow - vertdist/2;
  for j=1:Ncol
    xbeg    = sum(position{i}(1:(j-1))) + horidist/2;
    xend    = sum(position{i}(1:(j  ))) - horidist/2;
    pos(1) = xbeg*width  + horipos;
    pos(2) = ybeg*height + vertpos;
    pos(3) = (xend-xbeg)*width;
    pos(4) = (yend-ybeg)*height;
    h{i}{j} = uicontrol(fig, ...
      options{:}, ...
      'position', pos, ...
      'style',    style{i}{j}, ...
      'string',   string{i}{j}, ...
      'tag',      tag{i}{j}, ...
      'value',    value{i}{j}, ...
      'callback', callback{i}{j} ...
      );
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_creategui(hObject, eventdata, handles);
% define the position of each GUI element
position = {
  [2 1 1 1]
  [2 1 1 1]
  [2 1 1 1]
  [1]
  [1]
  [1]
  [1]
  [1 1]
  };

% define the style of each GUI element
style = {
  {'text' 'edit' 'edit' 'edit'}
  {'text' 'edit' 'edit' 'edit'}
  {'text' 'edit' 'edit' 'edit'}
  {'pushbutton'}
  {'pushbutton'}
  {'toggle'}
  {'toggle'}
  {'text' 'edit'}
  };

% define the descriptive string of each GUI element
string = {
  {'rotate'    0 0 0}
  {'translate' 0 0 0}
  {'scale'     1 1 1}
  {'redisplay'}
  {'apply'}
  {'toggle grid'}
  {'toggle axes'}
  {'alpha' 0.7}
  };

% define the value of each GUI element
value = {
  {[] [] [] []}
  {[] [] [] []}
  {[] [] [] []}
  {[]}
  {[]}
  {0}
  {0}
  {[] []}
  };

% define a tag for each GUI element
tag = {
  {'' 'rx' 'ry' 'rz'}
  {'' 'tx' 'ty' 'tz'}
  {'' 'sx' 'sy' 'sz'}
  {''}
  {''}
  {'toggle grid'}
  {'toggle axes'}
  {'' 'alpha'}
  };

% define the callback function of each GUI element
callback = {
  {[] @cb_redraw @cb_redraw @cb_redraw}
  {[] @cb_redraw @cb_redraw @cb_redraw}
  {[] @cb_redraw @cb_redraw @cb_redraw}
  {@cb_redraw}
  {@cb_apply}
  {@cb_redraw}
  {@cb_redraw}
  {[] @cb_redraw}
  };

fig = get(hObject, 'parent');
layoutgui(fig, [0.7 0.05 0.25 0.50], position, style, string, value, tag, callback);

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
  hs = triplot(template.elec.chanpos, [], [], 'nodes');
  set(hs, 'MarkerSize', 10);
  set(hs, 'Color', 'b');
  if isfield(template.elec, 'line')
    hs = triplot(template.elec.chanpos, template.elec.line, [], 'edges');
    try, set(hs, 'MarkerEdgeColor', 'b'); end
  end
end

if ~isempty(individual.elec)
  hs = triplot(individual.elec.chanpos, [], [], 'nodes');
  set(hs, 'MarkerSize', 10);
  set(hs, 'Color', 'r');
  if isfield(individual.elec, 'line')
    hs = triplot(individual.elec.chanpos, elec.line, [], 'edges');
    try, set(hs, 'MarkerEdgeColor', 'r'); end
  end
end

if ~isempty(template.grad)
  hs = triplot(template.grad.chanpos, [], [], 'nodes');
  set(hs, 'MarkerSize', 10);
  set(hs, 'Color', 'b');
  % FIXME also plot lines?
end

if ~isempty(individual.grad)
  hs = triplot(individual.grad.chanpos, [], [], 'nodes');
  set(hs, 'MarkerSize', 10);
  set(hs, 'Color', 'r');
  % FIXME also plot lines?
end

if ~isempty(template.vol)
  % FIXME this only works for boundary element models
  for i = 1:numel(template.vol.bnd)
    if strcmp(template.volstyle, 'edge') || ...
        strcmp(template.volstyle, 'both')
      hs = triplot(template.vol.bnd(i).pnt, template.vol.bnd(i).tri, [], 'edges');
      try, set(hs, 'EdgeColor', 'b'); end
    end
    if strcmp(template.volstyle, 'surface') || ...
        strcmp(template.volstyle, 'both')
      hs = triplot(template.vol.bnd(i).pnt, template.vol.bnd(i).tri, [], 'faces_blue');
      
    end
  end
end

if ~isempty(individual.vol)
  % FIXME this only works for boundary element models
  for i = 1:numel(individual.vol.bnd)
    if strcmp(individual.volstyle, 'edge') || ...
        strcmp(individual.volstyle, 'both')
      hs = triplot(individual.vol.bnd(i).pnt, individual.vol.bnd(i).tri, [], 'edges');
      try, set(hs, 'EdgeColor', 'r'); end
    end
    if strcmp(individual.volstyle, 'surface') || ...
        strcmp(individual.volstyle, 'both')
      hs = triplot(individual.vol.bnd(i).pnt, individual.vol.bnd(i).tri, [], 'faces_red');
    end
  end
end

if ~isempty(template.headshape)
  if isfield(template.headshape, 'pnt') && ~isempty(template.headshape.pnt)
    if strcmp(template.headshapestyle, 'surface') || ...
        strcmp(template.headshapestyle, 'both')
      triplot(template.headshape.pnt, template.headshape.tri,  [], 'faces_blue');
      alpha(str2num(get(findobj(fig, 'tag', 'alpha'), 'string')));
    end
    
    if strcmp(template.headshapestyle, 'vertex') || ...
        strcmp(template.headshapestyle, 'both')
      hs = triplot(template.headshape.pnt, [], [], 'nodes');
      set(hs, 'Color', 'b');
    end
  end
  if isfield(template.headshape, 'fid') && ~isempty(template.headshape.fid.pnt)
    triplot(template.headshape.fid.pnt, [], [], 'nodes_blue');
  end
end

if ~isempty(individual.headshape)
  if isfield(individual.headshape, 'pnt') && ~isempty(individual.headshape.pnt)
    if strcmp(individual.headshapestyle, 'surface') || ...
        strcmp(individual.headshapestyle, 'both')
      triplot(individual.headshape.pnt, individual.headshape.tri,  [], 'faces_red');
      alpha(str2num(get(findobj(fig, 'tag', 'alpha'), 'string')));
    end
    
    if strcmp(individual.headshapestyle, 'vertex') || ...
        strcmp(individual.headshapestyle, 'both')
      hs = triplot(individual.headshape.pnt, [], [], 'nodes');
      set(hs, 'Color', 'r');
    end
  end
  if isfield(individual.headshape, 'fid') && ~isempty(individual.headshape.fid.pnt)
    triplot(individual.headshape.fid.pnt, [], [], 'nodes_red');
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

