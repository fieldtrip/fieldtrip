function mri = ft_defacevolume(cfg, mri)

% FT_DEFACEVOLUME allows you to blank out specific regions from an anatomical MRI,
% such as the face and ears. The graphical user interface allows you to position a
% box over the anatomical MRI inside which all anatomical voxel values will be
% replaced by zero. Depending on the alignment of the anatomical MRI and whether both
% face and ears need to be removed, you might have to call this function multiple
% times in succession. Following defacing, you should check the result with
% FT_SOURCEPLOT.
%
% Use as
%   mri = ft_defacevolume(cfg, mri)
%
% The configuration can contain the following options
%   cfg.center = position of the initial box (default = [0 0 0])
%   cfg.size   = size of the initial box (default is automatic)
%
% See also FT_ANONIMIZEDATA, FT_ANALYSISPIPELINE, FT_SOURCEPLOT

% Copyright (C) 2015, Robert Oostenveld
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
ft_preamble loadvar    mri
ft_preamble provenance mri
ft_preamble trackconfig
ft_preamble debug

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% set the defaults
cfg.center = ft_getopt(cfg, 'center', [0 0 0]);
cfg.size   = ft_getopt(cfg, 'size'); % automatic default is determined further down

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');

% determine the size of the "unit" sphere in the origin and the length of the axes
switch mri.unit
  case 'mm'
    axmax = 150;
    rbol  = 5;
  case 'cm'
    axmax = 15;
    rbol  = 0.5;
  case 'm'
    axmax = 0.15;
    rbol  = 0.005;
  otherwise
    error('unknown units (%s)', unit);
end

% the volumetric data needs to be interpolated onto three orthogonal planes
% determine a resolution that is close to, or identical to the original resolution
[corner_vox, corner_head] = cornerpoints(mri.dim, mri.transform);
diagonal_head = norm(range(corner_head));
diagonal_vox  = norm(range(corner_vox));
resolution    = diagonal_head/diagonal_vox; % this is in units of "mri.unit"

figHandle = figure;
set(figHandle, 'CloseRequestFcn', @cb_close);

% clear persistent variables to ensure fresh figure
clear ft_plot_slice

% create a contrast enhanced version of the anatomy
mri.anatomy = double(mri.anatomy);
dum = unique(mri.anatomy(:));
clim(1) = dum(round(0.05*numel(dum)));
clim(2) = dum(round(0.95*numel(dum)));
anatomy = (mri.anatomy-clim(1))/(clim(2)-clim(1));

ft_plot_ortho(anatomy, 'transform', mri.transform, 'resolution', resolution, 'style', 'intersect');
axis vis3d
view([110 36]);

% shift the axes to the left
ax = get(gca, 'position');
ax(1) = 0;
set(gca, 'position', ax);

% get the xyz-axes
xdat  = [-axmax 0 0; axmax 0 0];
ydat  = [0 -axmax 0; 0 axmax 0];
zdat  = [0 0 -axmax; 0 0 axmax];

% get the xyz-axes dotted
xdatdot = (-axmax:(axmax/15):axmax);
xdatdot = xdatdot(1:floor(numel(xdatdot)/2)*2);
xdatdot = reshape(xdatdot, [2 numel(xdatdot)/2]);
n       = size(xdatdot,2);
ydatdot = [zeros(2,n) xdatdot zeros(2,n)];
zdatdot = [zeros(2,2*n) xdatdot];
xdatdot = [xdatdot zeros(2,2*n)];

% plot axes
hl = line(xdat, ydat, zdat);
set(hl(1), 'linewidth', 1, 'color', 'r');
set(hl(2), 'linewidth', 1, 'color', 'g');
set(hl(3), 'linewidth', 1, 'color', 'b');
hld = line(xdatdot, ydatdot, zdatdot);
for k = 1:n
  set(hld(k    ), 'linewidth', 3, 'color', 'r');
  set(hld(k+n*1), 'linewidth', 3, 'color', 'g');
  set(hld(k+n*2), 'linewidth', 3, 'color', 'b');
end

mask.pos = cfg.center;
mask.siz = cfg.size;
if isempty(mask.siz)
  mask.siz = [axmax axmax axmax]/2;
end

guidata(figHandle, mask);

% add the GUI elements
cb_creategui(gca);
cb_redraw(gca);

uiwait(figHandle);
mask = guidata(figHandle);
delete(figHandle);
drawnow
fprintf('applying mask to MRI voxels\n')

x1 = mask.pos(1) - mask.siz(1);
y1 = mask.pos(1) - mask.siz(1);
z1 = mask.pos(2) - mask.siz(2);
x2 = mask.pos(2) + mask.siz(2);
y2 = mask.pos(3) + mask.siz(3);
z2 = mask.pos(3) + mask.siz(3);

R = mask.R;
S = mask.S;
T = mask.T;

% it is possible to convert the box to headcoordinates, but it is more efficient the other way around
[X, Y, Z] = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));
voxpos = ft_warp_apply(mri.transform, [X(:) Y(:) Z(:)]);  % voxel positions in head coordinates
voxpos = ft_warp_apply(inv(T*S*R), voxpos);           % voxel positions in box coordinates

inside = ...
  voxpos(:,1) > -mask.siz(1) & ...
  voxpos(:,1) < +mask.siz(1) & ...
  voxpos(:,2) > -mask.siz(2) & ...
  voxpos(:,2) < +mask.siz(2) & ...
  voxpos(:,3) > -mask.siz(3) & ...
  voxpos(:,3) < +mask.siz(3);

mri.anatomy(inside) = 0;
fprintf('masked %.0f%% of the voxels\n', 100*mean(inside));

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble history mri
ft_postamble savevar mri


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(figHandle, varargin)
persistent p
% define the position of each GUI element
figHandle = get(figHandle, 'parent');
mask = guidata(figHandle);

rx = str2double(get(findobj(figHandle, 'tag', 'rx'), 'string'));
ry = str2double(get(findobj(figHandle, 'tag', 'ry'), 'string'));
rz = str2double(get(findobj(figHandle, 'tag', 'rz'), 'string'));
tx = str2double(get(findobj(figHandle, 'tag', 'tx'), 'string'));
ty = str2double(get(findobj(figHandle, 'tag', 'ty'), 'string'));
tz = str2double(get(findobj(figHandle, 'tag', 'tz'), 'string'));
sx = str2double(get(findobj(figHandle, 'tag', 'sx'), 'string'));
sy = str2double(get(findobj(figHandle, 'tag', 'sy'), 'string'));
sz = str2double(get(findobj(figHandle, 'tag', 'sz'), 'string'));

R = rotate   ([rx ry rz]);
T = translate([tx ty tz]);
S = scale    ([sx sy sz]);

x1 = mask.pos(1) - mask.siz(1);
y1 = mask.pos(1) - mask.siz(1);
z1 = mask.pos(2) - mask.siz(2);
x2 = mask.pos(2) + mask.siz(2);
y2 = mask.pos(3) + mask.siz(3);
z2 = mask.pos(3) + mask.siz(3);

plane1 = [
  x1 y1 z1
  x2 y1 z1
  x2 y2 z1
  x1 y2 z1];

plane2 = [
  x1 y1 z2
  x2 y1 z2
  x2 y2 z2
  x1 y2 z2];

plane3 = [
  x1 y1 z1
  x1 y2 z1
  x1 y2 z2
  x1 y1 z2];

plane4 = [
  x2 y1 z1
  x2 y2 z1
  x2 y2 z2
  x2 y1 z2];

plane5 = [
  x1 y1 z1
  x2 y1 z1
  x2 y1 z2
  x1 y1 z2];

plane6 = [
  x1 y2 z1
  x2 y2 z1
  x2 y2 z2
  x1 y2 z2];

plane1 = ft_warp_apply(T, ft_warp_apply(R*S, plane1));
plane2 = ft_warp_apply(T, ft_warp_apply(R*S, plane2));
plane3 = ft_warp_apply(T, ft_warp_apply(R*S, plane3));
plane4 = ft_warp_apply(T, ft_warp_apply(R*S, plane4));
plane5 = ft_warp_apply(T, ft_warp_apply(R*S, plane5));
plane6 = ft_warp_apply(T, ft_warp_apply(R*S, plane6));

if all(ishandle(p))
  delete(p);
end

p(1) = patch(plane1(:,1), plane1(:,2), plane1(:,3), 'y');
p(2) = patch(plane2(:,1), plane2(:,2), plane2(:,3), 'y');
p(3) = patch(plane3(:,1), plane3(:,2), plane3(:,3), 'y');
p(4) = patch(plane4(:,1), plane4(:,2), plane4(:,3), 'y');
p(5) = patch(plane5(:,1), plane5(:,2), plane5(:,3), 'y');
p(6) = patch(plane6(:,1), plane6(:,2), plane6(:,3), 'y');
set(p, 'FaceAlpha', 0.3);

mask.R = R;
mask.T = T;
mask.S = S;
guidata(figHandle, mask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_creategui(figHandle, varargin)
% define the position of each GUI element
figHandle = get(figHandle, 'parent');
% constants
CONTROL_WIDTH   = 0.05;
CONTROL_HEIGHT  = 0.08;
CONTROL_HOFFSET = 0.68;
CONTROL_VOFFSET = 0.20;

% rotateui
uicontrol('tag', 'rotateui', 'parent', figHandle, 'units', 'normalized', 'style', 'text', 'string', 'rotate', 'callback', [])
uicontrol('tag', 'rx', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'ry', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'rz', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
ft_uilayout(figHandle, 'tag', 'rotateui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET);
ft_uilayout(figHandle, 'tag', 'rx',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);
ft_uilayout(figHandle, 'tag', 'ry',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);
ft_uilayout(figHandle, 'tag', 'rz',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);

% scaleui
uicontrol('tag', 'scaleui', 'parent', figHandle, 'units', 'normalized', 'style', 'text', 'string', 'scale', 'callback', [])
uicontrol('tag', 'sx', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
uicontrol('tag', 'sy', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
uicontrol('tag', 'sz', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', '1', 'callback', @cb_redraw)
ft_uilayout(figHandle, 'tag', 'scaleui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'sx',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'sy',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'sz',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);

% translateui
uicontrol('tag', 'translateui', 'parent', figHandle, 'units', 'normalized', 'style', 'text', 'string', 'translate', 'callback', [])
uicontrol('tag', 'tx', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'ty', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
uicontrol('tag', 'tz', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', '0', 'callback', @cb_redraw)
ft_uilayout(figHandle, 'tag', 'translateui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'tx',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'ty',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'tz',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_close(figHandle, varargin)
% the figure will be closed in the main function after collecting the guidata
uiresume;