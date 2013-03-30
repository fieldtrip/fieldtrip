function [h, T2] = ft_plot_slice(dat, varargin)

% FT_PLOT_SLICE cuts a 2-D slice from a 3-D volume and interpolates if needed
%
% Use as
%   ft_plot_slice(dat, ...)
%   ft_plot_ortho(dat, mask, ...)
% where dat and mask are equal-sized 3-D arrays.
%
% Additional options should be specified in key-value pairs and can be
%   'transform'    = 4x4 homogeneous transformation matrix specifying the mapping from
%                    voxel space to the coordinate system in which the data are plotted.
%   'location'     = 1x3 vector specifying a point on the plane which will be plotted
%                    the coordinates are expressed in the coordinate system in which the
%                    data will be plotted. location defines the origin of the plane
%   'orientation'  = 1x3 vector specifying the direction orthogonal through the plane
%                    which will be plotted (default = [0 0 1])
%   'resolution'   = number (default = 1)
%   'datmask'      = 3D-matrix with the same size as the data matrix, serving as opacitymap
%                    If the second input argument to the function contains a matrix, this
%                    will be used as the mask
%   'opacitylim'   = 1x2 vector specifying the limits for opacity masking
%   'interpmethod' = string specifying the method for the interpolation, see INTERPN (default = 'nearest')
%   'style'        = string, 'flat' or '3D'
%   'colormap'     = string, see COLORMAP
%   'colorlim'     = 1x2 vector specifying the min and max for the colorscale
%
% See also FT_PLOT_ORTHO, FT_PLOT_MONTAGE, FT_SOURCEPLOT

% undocumented
%   'intersectmesh'  = triangulated mesh through which the intersection of the plane will be plotted (e.g. cortical sheet)
%   'intersectcolor' = color for the intersection

% Copyrights (C) 2010, Jan-Mathijs Schoffelen
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

persistent previous_dim X Y Z

% parse first input argument(s). it is either
% (dat, varargin)
% (dat, msk, varargin)
% (dat, [], varargin)
if numel(varargin)>0 && (isempty(varargin{1}) || isnumeric(varargin{1}) || islogical(varargin{1}))
  M        = varargin{1};
  varargin = varargin(2:end);
end

% get the optional input arguments
transform    = ft_getopt(varargin, 'transform');
loc          = ft_getopt(varargin, 'location');
ori          = ft_getopt(varargin, 'orientation', [0 0 1]);
resolution   = ft_getopt(varargin, 'resolution', 1);
mask         = ft_getopt(varargin, 'datmask');
opacitylim   = ft_getopt(varargin, 'opacitylim');
interpmethod = ft_getopt(varargin, 'interpmethod', 'nearest');
cmap         = ft_getopt(varargin, 'colormap');
clim         = ft_getopt(varargin, 'colorlim');
doscale      = ft_getopt(varargin, 'doscale', true); % only scale when necessary (time consuming), i.e. when plotting as grayscale image & when the values are not between 0 and 1
h            = ft_getopt(varargin, 'surfhandle', []);

mesh                = ft_getopt(varargin, 'intersectmesh');
intersectcolor      = ft_getopt(varargin, 'intersectcolor', 'yrgbmyrgbm');
intersectlinewidth  = ft_getopt(varargin, 'intersectlinewidth', 2);
intersectlinestyle  = ft_getopt(varargin, 'intersectlinestyle');

% convert from yes/no/true/false/0/1 into a proper boolean
doscale = istrue(doscale);

if ~isa(dat, 'double')
  dat = cast(dat, 'double');
end

if exist('M', 'var') && isempty(mask)
  warning_once('using the mask from the input and not from the varargin list');
  mask = M; clear M;
end

% norm normalise the ori vector
ori = ori./sqrt(sum(ori.^2));

% dimensionality of the input data
dim = size(dat);
if isempty(previous_dim),
  previous_dim = [0 0 0];
end

% set the location if empty
if isempty(loc) && (isempty(transform) || all(all(transform-eye(4)==0)==1))
  loc = dim./2;
elseif isempty(loc)
  loc = [0 0 0];
end

% set the transformation matrix if empty
if isempty(transform)
  transform = eye(4);
end

% check whether the mesh is ok
dointersect = ~isempty(mesh);
if ~iscell(mesh)
  mesh = {mesh};
end

if dointersect
  for k = 1:numel(mesh)
    if isfield(mesh{k}, 'pos')
      % use pos instead of pnt
      mesh{k}.pnt = mesh{k}.pos;
      mesh{k} = rmfield(mesh{k}, 'pos');
    end
    if ~isfield(mesh{k}, 'pnt') || ~isfield(mesh{k}, 'tri')
      error('the triangulated mesh should be a structure with pnt and tri');
    end
  end
end

% check whether the mask is ok
domask = ~isempty(mask);
if domask
  if ~isequal(size(dat), size(mask))
    error('the mask data should have the same dimensions as the functional data');
  end
end

% determine whether interpolation is needed
dointerp = false;
dointerp = dointerp || sum(sum(transform-eye(4)))~=0;
dointerp = dointerp || ~all(round(loc)==loc);
dointerp = dointerp || sum(ori)~=1;
dointerp = dointerp || ~(resolution==round(resolution));

% determine the caller function and toggle dointerp to true, if
% ft_plot_slice has been called from ft_plot_montage
% this is necessary for the correct allocation of the persistent variables
st = dbstack;
if ~dointerp && numel(st)>1 && strcmp(st(2).name, 'ft_plot_montage'), dointerp = true; end

% determine the corner points of the volume in voxel and in plotting space
[corner_vox, corner_head] = cornerpoints(dim+1, transform);
  
if dointerp
  %--------cut a slice using interpn
  
  % get voxel indices
  if all(dim==previous_dim)
    % for speeding up the plotting on subsequent calls
    % use persistent variables X Y Z
  else
    [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
  end
    
  % define 'x' and 'y' axis in projection plane, the definition of x and y is more or less arbitrary
  [x, y] = projplane(ori); % z = ori
  
  % project the corner points onto the projection plane
  corner_proj = nan(size(corner_head));
  for i=1:8
    corner   = corner_head(i, :);
    corner   = corner - loc(:)';
    corner_x = dot(corner, x);
    corner_y = dot(corner, y);
    corner_z = 0;
    corner_proj(i, :) = [corner_x corner_y corner_z];
  end
  
  % determine a tight grid of points in the projection plane
  xplane = floor(min(corner_proj(:, 1))):resolution:ceil(max(corner_proj(:, 1)));
  yplane = floor(min(corner_proj(:, 2))):resolution:ceil(max(corner_proj(:, 2)));
  zplane = 0;
  
  [X2, Y2, Z2] = ndgrid(xplane, yplane, zplane); %2D cartesian grid of projection plane in plane voxels
  siz        = size(squeeze(X2));
  pos        = [X2(:) Y2(:) Z2(:)]; clear X2 Y2 Z2;
  
  if false
    % this is for debugging
    ft_plot_mesh(warp_apply(T2, pos))
    ft_plot_mesh(corner_head)
    axis on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
  end
  
  % get the transformation matrix from plotting space to voxel space
  % T1 = inv(transform);
  
  % get the transformation matrix to get the projection plane at the right location and orientation into plotting space
  T2  = [x(:) y(:) ori(:) loc(:); 0 0 0 1];
  
  % get the transformation matrix from projection plane to voxel space
  % M = T1*T2;
  M = transform\T2;
  
  % get the positions of the pixels of the desires plane in voxel space
  pos = warp_apply(M, pos);
  
  Xi              = reshape(pos(:, 1), siz);
  Yi              = reshape(pos(:, 2), siz);
  Zi              = reshape(pos(:, 3), siz);
  V               = interpn(X, Y, Z, dat, Xi, Yi, Zi, interpmethod);
  [V, Xi, Yi, Zi] = tight(V, Xi, Yi, Zi);
  siz             = size(Xi);
  
  if domask,
    Vmask = tight(interpn(X, Y, Z, mask, Xi, Yi, Zi, interpmethod));
  end
  
  if ~isempty(Xi)
    % now adjust the Xi, Yi and Zi, to allow for the surface object
    % convention, where the data point value is defined in the center of each
    % square, i.e. the number of elements in each of the dimensions of Xi, Yi
    % Zi should be 1 more than the functional data, and they should be displaced
    % by half a voxel distance
    dx2 = mean(diff(Xi,[],2),2); dx1 = mean(diff(Xi,[],1),1);
    dy2 = mean(diff(Yi,[],2),2); dy1 = mean(diff(Yi,[],1),1);
    dz2 = mean(diff(Zi,[],2),2); dz1 = mean(diff(Zi,[],1),1);
    
    Xi  = [Xi-0.5*dx2*ones(1,siz(2)) Xi(:,end)+0.5*dx2; Xi(end,:)+0.5*dx1 Xi(end,end)+0.5*(dx1(end)+dx2(end))];
    Yi  = [Yi-0.5*dy2*ones(1,siz(2)) Yi(:,end)+0.5*dy2; Yi(end,:)+0.5*dy1 Yi(end,end)+0.5*(dy1(end)+dy2(end))];
    Zi  = [Zi-0.5*dz2*ones(1,siz(2)) Zi(:,end)+0.5*dz2; Zi(end,:)+0.5*dz1 Zi(end,end)+0.5*(dz1(end)+dz2(end))];
  end
  
else
  %-------cut a slice without interpolation
  [x, y] = projplane(ori);
  T2     = [x(:) y(:) ori(:) loc(:); 0 0 0 1];
  
  % the '+1' and '-0.5' are needed due to the difference between handling
  % of image and surf. Surf color data is defined in the center of each
  % square, hence needs axes and coordinate adjustment
  if all(ori==[1 0 0]), xplane = loc(1);   yplane = 1:(dim(2)+1); zplane = 1:(dim(3)+1); end
  if all(ori==[0 1 0]), xplane = 1:(dim(1)+1); yplane = loc(2);   zplane = 1:(dim(3)+1); end
  if all(ori==[0 0 1]), xplane = 1:(dim(1)+1); yplane = 1:(dim(2)+1); zplane = loc(3);   end
  
  [Xi, Yi, Zi] = ndgrid(xplane-0.5, yplane-0.5, zplane-0.5); % coordinate is centre of the voxel, 1-based
  siz        = size(squeeze(Xi))-1;
  Xi         = reshape(Xi, siz+1); if numel(xplane)==1, xplane = xplane([1 1]); end;
  Yi         = reshape(Yi, siz+1); if numel(yplane)==1, yplane = yplane([1 1]); end;
  Zi         = reshape(Zi, siz+1); if numel(zplane)==1, zplane = zplane([1 1]); end;
  V          = reshape(dat(xplane(1:end-1), yplane(1:end-1), zplane(1:end-1)), siz);
  if domask,
    Vmask    = reshape(mask(xplane(1:end-1), yplane(1:end-1), zplane(1:end-1)), siz);
  end
  
end

if isempty(cmap),
  % treat as gray value: scale and convert to rgb
  if doscale
    dmin = min(dat(:));
    dmax = max(dat(:));
    V    = (V-dmin)./(dmax-dmin);
    clear dmin dmax
  end
  V(isnan(V)) = 0;
  % convert anatomy into RGB values
  V = cat(3, V, V, V);
end

if isempty(h),
  % get positions of the plane in plotting space
  posh = warp_apply(transform, [Xi(:) Yi(:) Zi(:)], 'homogeneous', 1e-8);
  if ~isempty(posh)
    Xh   = reshape(posh(:, 1), siz+1);
    Yh   = reshape(posh(:, 2), siz+1);
    Zh   = reshape(posh(:, 3), siz+1);
  else
    % emulate old behavior, that allowed empty data to be plotted
    Xh   = [];
    Yh   = [];
    Zh   = [];
  end
  % create surface object
  h    = surface(Xh, Yh, Zh, V);
else
  % update the colordata in the surface object
  set(h, 'Cdata', V);
end

set(h, 'linestyle', 'none');

if domask,
  if islogical(Vmask), Vmask = double(Vmask); end
  set(h, 'FaceAlpha', 'flat');
  set(h, 'AlphaDataMapping', 'scaled');
  set(h, 'AlphaData', Vmask);
  if ~isempty(opacitylim)
    alim(opacitylim)
  end
end

if dointersect
  % determine three points on the plane
  inplane = eye(3) - (eye(3) * ori') * ori;
  v1 = loc + inplane(1,:);
  v2 = loc + inplane(2,:);
  v3 = loc + inplane(3,:);
  
  for k = 1:numel(mesh)
    [xmesh, ymesh, zmesh] = intersect_plane(mesh{k}.pnt, mesh{k}.tri, v1, v2, v3);
    
    % draw each individual line segment of the intersection
    if ~isempty(xmesh), 
      p(k) = patch(xmesh', ymesh', zmesh', nan(1, size(xmesh,1)));
      if ~isempty(intersectcolor),     set(p(k), 'EdgeColor', intersectcolor(k)); end
      if ~isempty(intersectlinewidth), set(p(k), 'LineWidth', intersectlinewidth); end
      if ~isempty(intersectlinestyle), set(p(k), 'LineStyle', intersectlinestyle); end
  end
  end
end

if ~isempty(cmap)
  colormap(cmap);
end

if ~isempty(clim)
  caxis(clim);
end

% update the axes to ensure that the whole volume fits
ax = [min(corner_head) max(corner_head)];
axis(ax([1 4 2 5 3 6])-0.5); % reorder into [xmin xmax ymin ymaz zmin zmax]
st = dbstack;
if numel(st)>1,
  % ft_plot_slice has been called from another function
  % assume the axis settings to be handled there
else
  axis equal
  axis vis3d
end
  
% store for future reference
previous_dim  = dim;

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
function [V, Xi, Yi, Zi] = tight(V, Xi, Yi, Zi)
% cut off the nans at the edges
x = sum(~isfinite(V), 1)<size(V, 1);
y = sum(~isfinite(V), 2)<size(V, 2);
V = V(y, x);
if nargin>1, Xi = Xi(y, x); end
if nargin>2, Yi = Yi(y, x); end
if nargin>3, Zi = Zi(y, x); end
