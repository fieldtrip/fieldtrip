function [h, T2] = ft_plot_slice(dat, varargin)

% FT_PLOT_SLICE plots a 2-D cut through a 3-D volume and interpolates if needed
%
% Use as
%   ft_plot_slice(dat, ...)
% or
%   ft_plot_slice(dat, mask, ...)
% where dat and mask are equal-sized 3-D arrays.
%
% Additional options should be specified in key-value pairs and can be
%   'transform'    = 4x4 homogeneous transformation matrix specifying the mapping from
%                    voxel coordinates to the coordinate system in which the data are plotted.
%   'location'     = 1x3 vector specifying a point on the plane which will be plotted
%                    the coordinates are expressed in the coordinate system in which the
%                    data will be plotted. location defines the origin of the plane
%   'orientation'  = 1x3 vector specifying the direction orthogonal through the plane
%                    which will be plotted (default = [0 0 1])
%   'unit'         = string, can be 'm', 'cm' or 'mm (default is automatic)
%   'resolution'   = number (default = 1 mm)
%   'datmask'      = 3D-matrix with the same size as the data matrix, serving as opacitymap
%                    If the second input argument to the function contains a matrix, this
%                    will be used as the mask
%   'maskstyle'    = string, 'opacity' or 'colormix', defines the rendering
%   'background'   = needed when maskstyle is 'colormix', 3D-matrix with
%                    the same size as the data matrix, serving as
%                    grayscale image that provides the background
%   'opacitylim'   = 1x2 vector specifying the limits for opacity masking
%   'interpmethod' = string specifying the method for the interpolation, see INTERPN (default = 'nearest')
%   'style'        = string, 'flat' or '3D'
%   'colormap'     = string, see COLORMAP
%   'clim'         = 1x2 vector specifying the min and max for the colorscale
%
% You can plot the slices from the volume together with an intersection of the slices
% with a triangulated surface mesh (e.g. a cortical sheet) using
%   'intersectmesh'       = triangulated mesh, see FT_PREPARE_MESH
%   'intersectcolor'      = string, color specification
%   'intersectlinestyle'  = string, line specification 
%   'intersectlinewidth'  = number
%
% See also FT_PLOT_ORTHO, FT_PLOT_MONTAGE, FT_SOURCEPLOT

% Undocumented options
%   'plotmarker'     = Nx3 matrix with points to be plotted as markers, e.g. dipole positions
%   'markersize'
%   'markercolor'

% Copyrights (C) 2010-2014, Jan-Mathijs Schoffelen
% Copyrights (C) 2014-2016, Robert Oostenveld and Jan-Mathijs Schoffelen
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

persistent dim X Y Z

if isequal(dim, size(dat(:,:,:,1,1)))
  % reuse the persistent variables to speed up subsequent calls with the same input
else
  dim       = size(dat); 
  if numel(dim)<3
    dim(3) = 1; % add 1 to catch size(dat,3) is singleton
  end
  [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
end

if any(dim==1)
  ft_error('it is not possible to plot a volume with a dimensionality of 1 in one of its dimensions');
end

% parse first input argument(s). it is either
% (dat, varargin)
% (dat, msk, varargin)
% (dat, [], varargin)
if numel(varargin)>0 && (isempty(varargin{1}) || isnumeric(varargin{1}) || islogical(varargin{1}))
  msk      = varargin{1};
  varargin = varargin(2:end);
end

% get the optional input arguments
transform           = ft_getopt(varargin, 'transform', eye(4));
loc                 = ft_getopt(varargin, 'location');
ori                 = ft_getopt(varargin, 'orientation', [0 0 1]);
unit                = ft_getopt(varargin, 'unit');       % the default will be determined further down
resolution          = ft_getopt(varargin, 'resolution'); % the default depends on the units and will be determined further down
datmask             = ft_getopt(varargin, 'datmask');
maskstyle           = ft_getopt(varargin, 'maskstyle', 'opacity');
background          = ft_getopt(varargin, 'background');
opacitylim          = ft_getopt(varargin, 'opacitylim');
interpmethod        = ft_getopt(varargin, 'interpmethod', 'nearest');
cmap                = ft_getopt(varargin, 'colormap');
clim                = ft_getopt(varargin, 'clim');
doscale             = ft_getopt(varargin, 'doscale', true); % only scale when necessary (time consuming), i.e. when plotting as grayscale image & when the values are not between 0 and 1
h                   = ft_getopt(varargin, 'surfhandle', []);
p                   = ft_getopt(varargin, 'patchhandle', []);

mesh                = ft_getopt(varargin, 'intersectmesh');
intersectcolor      = ft_getopt(varargin, 'intersectcolor', 'yrgbmyrgbm');
intersectlinewidth  = ft_getopt(varargin, 'intersectlinewidth', 2);
intersectlinestyle  = ft_getopt(varargin, 'intersectlinestyle');

plotmarker          = ft_getopt(varargin, 'plotmarker');
markersize          = ft_getopt(varargin, 'markersize', 'auto');
markercolor         = ft_getopt(varargin, 'markercolor', 'w');

% convert from yes/no/true/false/0/1 into a proper boolean
doscale = istrue(doscale);

if ~isa(dat, 'double')
  dat = cast(dat, 'double');
end

if exist('msk', 'var') && isempty(datmask)
  ft_warning('using the second input argument as mask rather than the one from the varargin list');
  datmask = msk; clear msk;
end

% normalise the orientation vector to one
ori = ori./sqrt(sum(ori.^2));

% set the default location
if isempty(loc) && (isempty(transform) || isequal(transform, eye(4)))
  loc = (dim+1)./2;
elseif isempty(loc)
  loc = [0 0 0];
end

% shift the location to be along the orientation vector
loc = ori*dot(loc,ori);

% it should be a cell-array
if isstruct(mesh)
  tmp = mesh;
  mesh = cell(size(tmp));
  for i=1:numel(tmp)
    mesh{i} = tmp(i);
  end
elseif iscell(mesh)
  % do nothing
else
  mesh = {};
end

% replace pnt by pos
for k = 1:numel(mesh)
  mesh{k} = fixpos(mesh{k});
end

dointersect = ~isempty(mesh);
if dointersect
  for k = 1:numel(mesh)
    if ~isfield(mesh{k}, 'pos') || ~isfield(mesh{k}, 'tri')
      % ft_error('the mesh should be a structure with pos and tri');
      mesh{k}.pos = [];
      mesh{k}.tri = [];
    end
  end
end

% check whether the mask is ok
domask = ~isempty(datmask);
if domask
  if ~isequal(size(dat), size(datmask)) && ~isequal(cmap, 'rgb')
    % the exception is when the functional data is to be interpreted as rgb
    ft_error('the mask data should have the same dimensions as the functional data');
  end
end

dobackground = ~isempty(background);
if dobackground
  if ~isequal(size(dat), size(background))
    error('the background data should have the same dimensions as the functional data');
  end
end

% determine the voxel center
% voxel_center_vc = [X(:) Y(:) Z(:)];
% voxel_center_hc = ft_warp_apply(transform, voxel_center_vc);

% determine the edges, i.e. the corner points of each voxel
% [Xe, Ye, Ze] = ndgrid(0:dim(1), 0:dim(2), 0:dim(3));
% Xe = Xe+0.5;
% Ye = Ye+0.5;
% Ze = Ze+0.5;
% voxel_edge_vc = [Xe(:) Ye(:) Ze(:)];
% voxel_edge_hc = ft_warp_apply(transform, voxel_edge_vc);

% determine the corner points of the box encompassing the whole data block
% extend the box with half a voxel  in all directions to get the outer edge
corner_vc = [
  0.5        0.5        0.5
  0.5+dim(1) 0.5        0.5
  0.5+dim(1) 0.5+dim(2) 0.5
  0.5        0.5+dim(2) 0.5
  0.5        0.5        0.5+dim(3)
  0.5+dim(1) 0.5        0.5+dim(3)
  0.5+dim(1) 0.5+dim(2) 0.5+dim(3)
  0.5        0.5+dim(2) 0.5+dim(3)
  ];
corner_hc = ft_warp_apply(transform, corner_vc);

if isempty(unit)
  if ~isequal(transform, eye(4))
    % estimate the geometrical units we are dealing with
    unit = ft_estimate_units(norm(range(corner_hc)));
  else
    % units are in voxels, these are assumed to be close to mm
    unit = 'mm';
  end
end
if isempty(resolution)
  % the default resolution is 1 mm
  resolution = ft_scalingfactor('mm', unit);
end

% determine whether interpolation is needed
dointerp = false;
dointerp = dointerp || sum(sum(transform-eye(4)))~=0;
dointerp = dointerp || ~all(round(loc)==loc);
dointerp = dointerp || sum(ori)~=1;
dointerp = dointerp || ~(resolution==round(resolution));
% determine the caller function and toggle dointerp to true, if ft_plot_slice has been called from ft_plot_montage
% this is necessary for the correct allocation of the persistent variables
st = dbstack;
if ~dointerp && numel(st)>1 && strcmp(st(2).name, 'ft_plot_montage'), dointerp = true; end



% define 'x' and 'y' axis in projection plane, the definition of x and y is more or less arbitrary
[x, y] = projplane(ori);
% z = ori;

% project the corner points onto the projection plane
corner_pc = zeros(size(corner_hc));
for i=1:8
  corner   = corner_hc(i, :) - loc(:)';
  corner_pc(i,1) = dot(corner, x);
  corner_pc(i,2) = dot(corner, y);
  corner_pc(i,3) = 0;
end

% get the transformation matrix from the projection plane to head coordinates
T2 = [x(:) y(:) ori(:) loc(:); 0 0 0 1];

% get the transformation matrix from projection plane to voxel coordinates
T3 = transform\T2;

min_corner_pc = min(corner_pc, [], 1);
max_corner_pc = max(corner_pc, [], 1);
% round the bounding box limits to the nearest mm
switch unit
  case 'm'
    min_corner_pc = ceil(min_corner_pc*100)/100;
    max_corner_pc = floor(max_corner_pc*100)/100;
  case 'cm'
    min_corner_pc = ceil(min_corner_pc*10)/10;
    max_corner_pc = floor(max_corner_pc*10)/10;
  case 'mm'
    min_corner_pc = ceil(min_corner_pc);
    max_corner_pc = floor(max_corner_pc);
end

% determine a grid of points in the projection plane
xplane = min_corner_pc(1):resolution:max_corner_pc(1);
yplane = min_corner_pc(2):resolution:max_corner_pc(2);
zplane = 0;
[Xi, Yi, Zi]      = ndgrid(xplane, yplane, zplane);
siz               = [size(squeeze(Xi)) size(dat,4)];
interp_center_pc  = [Xi(:) Yi(:) Zi(:)];
% interp_center_hc = ft_warp_apply(T2, interp_center_pc);

% get the positions of the points in the projection plane in voxel coordinates
interp_center_vc = ft_warp_apply(T3, interp_center_pc);

Xi = reshape(interp_center_vc(:, 1), siz(1:2));
Yi = reshape(interp_center_vc(:, 2), siz(1:2));
Zi = reshape(interp_center_vc(:, 3), siz(1:2));

% check whether the values in the axes are close enough to integer
tol = nanmean([diff(unique(Xi(:)));diff(unique(Yi(:)))])./100;
isintegerXi = issufficientlyinteger(Xi(:),tol);
isintegerYi = issufficientlyinteger(Yi(:),tol);
isintegerZi = issufficientlyinteger(Zi(:),tol);

% check whether it's possible to select an orthogonal plane
[islineXi, lineXi] = isline(Xi);
[islineYi, lineYi] = isline(Yi);
[islineZi, lineZi] = isline(Zi);

use_interpn = ~isequal(transform, eye(4)) || ~isequal(interpmethod, 'nearest') || ~all([isintegerXi isintegerYi isintegerZi]);
get_slice   = ~use_interpn && all([islineXi islineYi islineZi]) && all([isintegerXi isintegerYi isintegerZi]);
if use_interpn
  V  = interpn(X, Y, Z, dat, Xi, Yi, Zi, interpmethod);
  if domask,       Vmask = interpn(X, Y, Z, datmask,    Xi, Yi, Zi, interpmethod); end
  if dobackground, Vback = interpn(X, Y, Z, background, Xi, Yi, Zi, interpmethod); end
elseif get_slice 
  %something more efficient than an interpolation can be done
  % just select the appropriate plane, and permute to get the orientation
  % right in the plots, something to do with ndgrid vs meshgrid I think
  permutevec = [2 1];
  if ndims(dat)>3
    permutevec = [permutevec 3:ndims(dat)];
  end
  if numel(unique(lineXi(:)))==1
    lineXi = lineXi(1);
  elseif numel(unique(lineYi(:)))==1
    lineYi = lineYi(1);
  elseif numel(unique(lineZi(:)))==1
    lineZi = lineZi(1);
  end
  V = permute(reshape(dat(lineXi,lineYi,lineZi,:), siz(permutevec(1:ndims(dat)-1))), permutevec);
  if domask,       Vmask = permute(reshape(datmask(lineXi,lineYi,lineZi,:),       siz(permutevec(1:2))), [2 1]); end
  if dobackground, Vback = permute(reshape(background(lineXi,lineYi,lineZi,:), siz(permutevec(1:2))), [2 1]); end
else
  % use sub2ind in the unlikely case that it's an oblique plane, parallel
  % to one of the axes with only integer indices
  % this fails for rgb data
  V = dat(sub2ind(dim, Xi(:), Yi(:), Zi(:)));
  V = reshape(V, siz);
end

if all(isnan(V(:)))
  % the projection plane lies completely outside the box spanned by the data
else
  % trim the edges of the projection plane
  [sel1, sel2] = tight(V(:,:,1));
  V  = V (sel1,sel2,:);
  Xi = Xi(sel1,sel2);
  Yi = Yi(sel1,sel2);
  Zi = Zi(sel1,sel2);
  if domask
    Vmask = Vmask(sel1,sel2);
  end
  if dobackground
    Vback = Vback(sel1,sel2);
  end
end

if dobackground
  % convert the background plane to a grayscale image
  bmin  = nanmin(background(:));
  bmax  = nanmax(background(:));
  Vback = (Vback-bmin)./(bmax-bmin);
  Vback(~isfinite(Vback)) = 0;
  Vback = cat(3, Vback, Vback, Vback);
end

interp_center_vc = [Xi(:) Yi(:) Zi(:)]; clear Xi Yi Zi
interp_center_pc = ft_warp_apply(inv(T3), interp_center_vc);

% determine a grid of points in the projection plane
% this reconstruction is needed since the edges may have been trimmed off
xplane = min(interp_center_pc(:, 1)):resolution:max(interp_center_pc(:, 1));
yplane = min(interp_center_pc(:, 2)):resolution:max(interp_center_pc(:, 2));
zplane = 0;

[Xi, Yi, Zi] = ndgrid(xplane, yplane, zplane); % 2D cartesian grid of projection plane in plane voxels
siz          = size(squeeze(Xi));

% extend with one voxel along dim 1
Xi = cat(1, Xi, Xi(end,:)+mean(diff(Xi,[],1),1));
Yi = cat(1, Yi, Yi(end,:)+mean(diff(Yi,[],1),1));
Zi = cat(1, Zi, Zi(end,:)+mean(diff(Zi,[],1),1));
% extend with one voxel along dim 2
Xi = cat(2, Xi, Xi(:,end)+mean(diff(Xi,[],2),2));
Yi = cat(2, Yi, Yi(:,end)+mean(diff(Yi,[],2),2));
Zi = cat(2, Zi, Zi(:,end)+mean(diff(Zi,[],2),2));
% shift with half a voxel along dim 1 and 2
Xi = Xi-0.5*resolution;
Yi = Yi-0.5*resolution;
% Zi = Zi; % do not shift along this direction

interp_edge_pc = [Xi(:) Yi(:) Zi(:)]; clear Xi Yi Zi
interp_edge_hc = ft_warp_apply(T2, interp_edge_pc);

if false
  % plot all objects in head coordinates
  ft_plot_mesh(voxel_center_hc,   'vertexmarker', 'o')
  ft_plot_mesh(voxel_edge_hc,     'vertexmarker', '+')
  ft_plot_mesh(corner_hc,         'vertexmarker', '*')
  ft_plot_mesh(interp_center_hc,  'vertexmarker', 'o', 'vertexcolor', 'r')
  ft_plot_mesh(interp_edge_hc,    'vertexmarker', '+', 'vertexcolor', 'r')
  axis on
  grid on
  xlabel('x')
  ylabel('y')
  zlabel('z')
end

if isempty(cmap)
  % treat as gray value: scale and convert to rgb
  if doscale
    dmin = min(dat(:));
    dmax = max(dat(:));
    V    = (V-dmin)./(dmax-dmin);
    clear dmin dmax
  end
  V(~isfinite(V)) = 0;
  
  % deal with clim for RGB data here, where the purpose is to increase the
  % contrast range, rather than shift the average grey value
  if ~isempty(clim)
    V = (V-clim(1))./clim(2);
    V(V>1)=1;
  end
  
  % convert into RGB values, e.g. for the plotting of anatomy
  V = cat(3, V, V, V);
end

% get positions of the voxels in the interpolation plane in head coordinates
Xh = reshape(interp_edge_hc(:,1), siz+1);
Yh = reshape(interp_edge_hc(:,2), siz+1);
Zh = reshape(interp_edge_hc(:,3), siz+1);

% do the actual plotting of the slice
if ~domask
  % no masked slice to be plotted
  if isempty(h)
    % create surface object
    h = surface(Xh, Yh, Zh, V);
    set(h, 'linestyle', 'none');
  else
    % update the colordata in the surface object
    set(h, 'Cdata', V);
    set(h, 'Xdata', Xh);
    set(h, 'Ydata', Yh);
    set(h, 'Zdata', Zh);
  end
elseif domask
  % what should be done depends on the maskstyle
  switch maskstyle
    case 'opacity'
      if dobackground
        ft_warning('specifying maskstyle = ''opacity'' causes the supplied background image not to be used');
      end
      if isempty(h)
        % create surface object
        h = surface(Xh, Yh, Zh, V);
        set(h, 'linestyle', 'none');
      else
        % update the colordata in the surface object
        set(h, 'Cdata', V);
        set(h, 'Xdata', Xh);
        set(h, 'Ydata', Yh);
        set(h, 'Zdata', Zh);
      end
      if islogical(Vmask), Vmask = double(Vmask); end
      set(h, 'FaceColor', 'texture');
      set(h, 'FaceAlpha', 'texturemap'); %flat
      set(h, 'AlphaDataMapping', 'scaled');
      set(h, 'AlphaData', Vmask);
      if ~isempty(opacitylim)
        alim(opacitylim)
      end
    
    case 'colormix'
      if isempty(cmap), error('using ''colormix'' as maskstyle requires an explicitly defined colormap'); end
      V = bg_rgba2rgb(Vback,V,cmap,clim,Vmask,'rampup',opacitylim);
      if isempty(h)
        % create surface object
        h = surface(Xh, Yh, Zh, V);
        set(h, 'linestyle', 'none');
      else
        % update the colordata in the surface object
        set(h, 'Cdata', V);
        set(h, 'Xdata', Xh);
        set(h, 'Ydata', Yh);
        set(h, 'Zdata', Zh);
      end
  otherwise
    error('unsupported maskstyle');
  end
end

% plot the intersection with a mesh
if dointersect
  % determine three points on the plane
  inplane = eye(3) - (eye(3) * ori') * ori;
  v1 = loc + inplane(1,:);
  v2 = loc + inplane(2,:);
  v3 = loc + inplane(3,:);
  
  for k = 1:numel(mesh)
    [xmesh, ymesh, zmesh] = intersect_plane(mesh{k}.pos, mesh{k}.tri, v1, v2, v3);
    
    % draw each individual line segment of the intersection
    if ~isempty(xmesh)
      if isempty(p)
        p = patch(xmesh', ymesh', zmesh', nan(1, size(xmesh,1)));
        if ~isempty(intersectcolor),     set(p, 'EdgeColor', intersectcolor(k));  end
        if ~isempty(intersectlinewidth), set(p, 'LineWidth', intersectlinewidth); end
        if ~isempty(intersectlinestyle), set(p, 'LineStyle', intersectlinestyle); end
      else
        set(p, 'XData', xmesh', 'YData', ymesh', 'ZData', zmesh', 'FaceVertexCdata', nan(size(xmesh,1),1));
      end
    end
  end
end

if ~isempty(cmap) && ~isequal(cmap, 'rgb')
  ft_colormap(cmap);
  if ~isempty(clim)
    caxis(clim);
  end
end

if ~isempty(plotmarker)
  % determine three points on the plane
  inplane = eye(3) - (eye(3) * ori') * ori;
  v1 = loc + inplane(1,:);
  v2 = loc + inplane(2,:);
  v3 = loc + inplane(3,:);
  pr = nan(size(plotmarker,1), 3);
  d  = nan(size(plotmarker,1), 1);
  for k = 1:size(plotmarker,1)
    [pr(k,:), d(k,:)] = ptriprojn(v1, v2, v3, plotmarker(k,:));
  end
  sel = d<eps*1e8;
  if sum(sel)>0
    ft_plot_dipole(pr(sel,:), repmat([0;0;1], 1, size(pr,1)), 'length', 0, 'color', markercolor, 'diameter', markersize);
  end
end

% update the axes to ensure that the whole volume fits
ax = [min(corner_hc) max(corner_hc)];
axis(ax([1 4 2 5 3 6])); % reorder into [xmin xmax ymin ymaz zmin zmax]

st = dbstack;
if numel(st)>1
  % ft_plot_slice has been called from another function
  % assume the remainder of the axis settings to be handled there
else
  set(gca,'xlim',[min(Xh(:))-0.5*resolution max(Xh(:))+0.5*resolution]);
  set(gca,'ylim',[min(Yh(:))-0.5*resolution max(Yh(:))+0.5*resolution]);
  set(gca,'zlim',[min(Zh(:))-0.5*resolution max(Zh(:))+0.5*resolution]);
  
  set(gca,'dataaspectratio',[1 1 1]);
  % axis equal; % this for some reason does not work robustly when drawing intersections, replaced by the above
  axis vis3d
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y] = projplane(z)
[u, s, v] = svd([eye(3) z(:)]);
x = u(:, 2)';
y = u(:, 3)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sel1, sel2] = tight(V)
% make a selection to cut off the nans at the edges
sel1 = sum(~isfinite(V), 2)<size(V, 2);
sel2 = sum(~isfinite(V), 1)<size(V, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = issufficientlyinteger(X, tolerance)
%isinteger only checks for integer class, so will always return false with
%double integers
bool = all(abs(X-round(X))<tolerance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bool, lineX] = isline(X)
%isline returns an array of the values in X are columnwise or rowwise the
%same, otherwise returns false
if isequal(X(ones(1,size(X,1)),:),X)
  lineX = X(1,:);
  bool = true;
elseif isequal(X(:,ones(1,size(X,2))),X)
  lineX = X(:,1)';
  bool = true;
else
  lineX = [];
  bool = false;
end

