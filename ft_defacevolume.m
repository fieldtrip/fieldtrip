function mri = ft_defacevolume(cfg, mri)

% FT_DEFACEVOLUME allows you to erase specific regions from an anatomical MRI,
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
%   cfg.translate  = initial position of the center of the box (default = [0 0 0])
%   cfg.scale      = initial size of the box along each dimension (default is automatic)
%   cfg.translate  = initial rotation of the box (default = [0 0 0])
%   cfg.selection  = which voxels to keep, can be 'inside' or 'outside' (default = 'outside')
%   cfg.smooth     = 'no' or the FWHM of the gaussian kernel in voxels (default = 'no')
%
% If you specify no smoothing, the selected area will be zero-masked. If you
% specify a certain amount of smoothing (in voxels FWHM), the selected area will
% be replaced by a smoothed version of the data.
%
% See also FT_ANONIMIZEDATA, FT_ANALYSISPIPELINE, FT_SOURCEPLOT

% Copyright (C) 2015, Robert Oostenveld
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
ft_preamble loadvar    mri
ft_preamble provenance mri
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.rotate    = ft_getopt(cfg, 'rotate', [0 0 0]);
cfg.scale     = ft_getopt(cfg, 'scale'); % the automatic default is determined further down
cfg.translate = ft_getopt(cfg, 'translate', [0 0 0]);
cfg.selection = ft_getopt(cfg, 'selection', 'outside');
cfg.smooth    = ft_getopt(cfg, 'smooth', 'no');

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

if isempty(cfg.scale)
  cfg.scale = [axmax axmax axmax]/2;
end

guidata(figHandle, cfg);

% add the GUI elements
cb_creategui(gca);
cb_redraw(gca);

uiwait(figHandle);
cfg = guidata(figHandle);
delete(figHandle);
drawnow
fprintf('keeping all voxels from MRI that are %s the box\n', cfg.selection)

R = cfg.R;
S = cfg.S;
T = cfg.T;

% it is possible to convert the box to headcoordinates, but it is more efficient the other way around
[X, Y, Z] = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));
voxpos = ft_warp_apply(mri.transform, [X(:) Y(:) Z(:)]);  % voxel positions in head coordinates
voxpos = ft_warp_apply(inv(T*S*R), voxpos);               % voxel positions in box coordinates

keep = ...
  voxpos(:,1) > -0.5 & ...
  voxpos(:,1) < +0.5 & ...
  voxpos(:,2) > -0.5 & ...
  voxpos(:,2) < +0.5 & ...
  voxpos(:,3) > -0.5 & ...
  voxpos(:,3) < +0.5;

if strcmp(cfg.selection, 'inside')
  % invert the selection, i.e. keep the voxels inside the box
  keep = ~keep;
end

if isequal(cfg.smooth, 'no')
  fprintf('zero-filling %.0f%% of the volume\n', 100*mean(keep));
  mri.anatomy(keep) = 0;
else
  tmp = mri.anatomy;
  tmp = (1 + 0.5.*randn(size(tmp))).*tmp; % add 50% noise to each voxel
  tmp = volumesmooth(tmp, cfg.smooth, 'anatomy');
  fprintf('smoothing %.0f%% of the volume\n', 100*mean(keep));
  mri.anatomy(keep) = tmp(keep);
end

% remove the temporary fields from the configuration, keep the rest for provenance
cfg = removefields(cfg, {'R', 'S', 'T'});

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous mri
ft_postamble provenance mri
ft_postamble history mri
ft_postamble savevar mri


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(figHandle, varargin)
persistent p
% define the position of each GUI element
figHandle = get(figHandle, 'parent');
cfg = guidata(figHandle);

rx = str2double(get(findobj(figHandle, 'tag', 'rx'), 'string'));
ry = str2double(get(findobj(figHandle, 'tag', 'ry'), 'string'));
rz = str2double(get(findobj(figHandle, 'tag', 'rz'), 'string'));
tx = str2double(get(findobj(figHandle, 'tag', 'tx'), 'string'));
ty = str2double(get(findobj(figHandle, 'tag', 'ty'), 'string'));
tz = str2double(get(findobj(figHandle, 'tag', 'tz'), 'string'));
sx = str2double(get(findobj(figHandle, 'tag', 'sx'), 'string'));
sy = str2double(get(findobj(figHandle, 'tag', 'sy'), 'string'));
sz = str2double(get(findobj(figHandle, 'tag', 'sz'), 'string'));

% remember the user specified transformation
cfg.rotate    = [rx ry rz];
cfg.translate = [tx ty tz];
cfg.scale     = [sx sy sz];

R = rotate   (cfg.rotate);
T = translate(cfg.translate);
S = scale    (cfg.scale);

% remember the transformation matrices
cfg.R = R;
cfg.T = T;
cfg.S = S;

% start with a cube of unit dimensions
x1 = -0.5;
y1 = -0.5;
z1 = -0.5;
x2 = +0.5;
y2 = +0.5;
z2 = +0.5;

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

plane1 = ft_warp_apply(T*R*S, plane1);
plane2 = ft_warp_apply(T*R*S, plane2);
plane3 = ft_warp_apply(T*R*S, plane3);
plane4 = ft_warp_apply(T*R*S, plane4);
plane5 = ft_warp_apply(T*R*S, plane5);
plane6 = ft_warp_apply(T*R*S, plane6);

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

guidata(figHandle, cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_creategui(figHandle, varargin)
% define the position of each GUI element
figHandle = get(figHandle, 'parent');
cfg = guidata(figHandle);

% constants
CONTROL_WIDTH   = 0.05;
CONTROL_HEIGHT  = 0.08;
CONTROL_HOFFSET = 0.68;
CONTROL_VOFFSET = 0.20;

% rotateui
uicontrol('tag', 'rotateui', 'parent', figHandle, 'units', 'normalized', 'style', 'text', 'string', 'rotate', 'callback', [])
uicontrol('tag', 'rx', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.rotate(1)), 'callback', @cb_redraw)
uicontrol('tag', 'ry', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.rotate(2)), 'callback', @cb_redraw)
uicontrol('tag', 'rz', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.rotate(3)), 'callback', @cb_redraw)
ft_uilayout(figHandle, 'tag', 'rotateui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET);
ft_uilayout(figHandle, 'tag', 'rx',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);
ft_uilayout(figHandle, 'tag', 'ry',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);
ft_uilayout(figHandle, 'tag', 'rz',       'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET);

% scaleui
uicontrol('tag', 'scaleui', 'parent', figHandle, 'units', 'normalized', 'style', 'text', 'string', 'scale', 'callback', [])
uicontrol('tag', 'sx', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.scale(1)), 'callback', @cb_redraw)
uicontrol('tag', 'sy', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.scale(2)), 'callback', @cb_redraw)
uicontrol('tag', 'sz', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.scale(3)), 'callback', @cb_redraw)
ft_uilayout(figHandle, 'tag', 'scaleui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'sx',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'sy',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'sz',      'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-CONTROL_HEIGHT);

% translateui
uicontrol('tag', 'translateui', 'parent', figHandle, 'units', 'normalized', 'style', 'text', 'string', 'translate', 'callback', [])
uicontrol('tag', 'tx', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.translate(1)), 'callback', @cb_redraw)
uicontrol('tag', 'ty', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.translate(2)), 'callback', @cb_redraw)
uicontrol('tag', 'tz', 'parent', figHandle, 'units', 'normalized', 'style', 'edit', 'string', num2str(cfg.translate(3)), 'callback', @cb_redraw)
ft_uilayout(figHandle, 'tag', 'translateui', 'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET,                 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'tx',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+3*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'ty',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+4*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);
ft_uilayout(figHandle, 'tag', 'tz',          'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH,   'height', CONTROL_HEIGHT/2, 'hpos', CONTROL_HOFFSET+5*CONTROL_WIDTH, 'vpos', CONTROL_VOFFSET-2*CONTROL_HEIGHT);

% somehow the toolbar gets lost in 2012b
set(figHandle, 'toolbar', 'figure');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_close(figHandle, varargin)
% the figure will be closed in the main function after collecting the guidata
uiresume;
