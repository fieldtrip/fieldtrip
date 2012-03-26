function [h, T2] = ft_plot_slice(dat, varargin)

% FT_PLOT_SLICE cuts a 2-D slice from a 3-D volume and interpolates if needed
%
% Use as
%   ft_plot_slice(dat, ...)
% or
%   ft_plot_ortho(dat, mask, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'transform'    = 4x4 homogeneous transformation matrix specifying the mapping from
%                    voxel space to the coordinate system in which the data are plotted.
%   'location'     = 1x3 vector specifying a point on the plane which will be plotted
%                    the coordinates are expressed in the coordinate system in which the
%                    data will be plotted. location defines the origin of the plane
%   'orientation'  = 1x3 vector specifying the direction orthogonal through the plane
%                    which will be plotted
%   'resolution'   = number (default = 1)
%   'datmask'      = 3D-matrix with the same size as the matrix dat, serving as opacitymap
%                    if the second input argument to the function
%                    contains a matrix, this will be used as the mask
%   'opacitylim'   = 1x2 vector specifying the limits for opacity masking
%   'interpmethod' = string specifying the method for the interpolation, 
%                    see INTERPN (default = 'nearest')
%   'style'        = string, 'flat' or '3D'
%   'colormap'     = string, see COLORMAP
%   'colorlim'     = 1x2 vector specifying the min and max for the
%                     colorscale
%   'interplim'
%
% See also FT_PLOT_ORTHO, FT_SOURCEPLOT

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

persistent previous_dim X Y Z;

% parse first input argument(s). it is either
% (dat, varargin)
% (dat, msk, varargin)
% (dat, [], varargin)
if isempty(varargin{1}) || isnumeric(varargin{1})
  M        = varargin{1};
  varargin = varargin(2:end);
end

% get the optional input arguments
transform    = ft_getopt(varargin, 'transform');
loc          = ft_getopt(varargin, 'location');
ori          = ft_getopt(varargin, 'orientation',  [0 0 1]);
resolution   = ft_getopt(varargin, 'resolution',   1);
mask         = ft_getopt(varargin, 'datmask');
opacitylim   = ft_getopt(varargin, 'opacitylim');
interpmethod = ft_getopt(varargin, 'interpmethod', 'nearest');
cmap         = ft_getopt(varargin, 'colormap');
clim         = ft_getopt(varargin, 'colorlim'); 
doscale      = ft_getopt(varargin, 'doscale', true); % only scale when necessary (time consuming), i.e. when plotting as grayscale image & when the values are not between 0 and 1
h            = ft_getopt(varargin, 'surfhandle', []);

doscale = istrue(doscale);

if ~isa(dat, 'double')
  dat = cast(dat, 'double');
end

if exist('M', 'var') && isempty(mask)
  warning_once('using the mask from the input and not from the varargin list');
  mask = M;clear M;
end

% norm normalise the ori vector
ori = ori./sqrt(sum(ori.^2));

% dimensionality of the input data
dim          = size(dat);
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

% check whether mask is ok
domask = ~isempty(mask);
if domask,
  if ~all(dim==size(mask)),
    error('the mask data should have the same dimensionality as the functional data');
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

if dointerp
  %--------cut a slice using interpn
  
  % get voxel indices
  if all(dim==previous_dim)
    % for speeding up the plotting on subsequent calls
    % use persistent variables X Y Z
  else
    [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
  end
  
  % determine the corner points of the volume in voxel and in plotting space
  [corner_vox, corner_head] = cornerpoints(dim, transform);
  
  % define 'x' and 'y' axis in projection plane, the definition of x and y is more or less arbitrary
  [x, y] = projplane(ori); % z = ori
  
  % project the corner points onto the projection plane
  corner_proj = nan(size(corner_head));
  for i=1:8
    corner = corner_head(i,:);
    corner = corner - loc(:)';
    corner_x = dot(corner, x);
    corner_y = dot(corner, y);
    corner_z = 0;
    corner_proj(i,:) = [corner_x corner_y corner_z];
  end
  
  % determine a tight grid of points in the projection plane
  xplane = floor(min(corner_proj(:,1))):resolution:ceil(max(corner_proj(:,1)));
  yplane = floor(min(corner_proj(:,2))):resolution:ceil(max(corner_proj(:,2)));
  zplane = 0;
  
  [X2,Y2,Z2] = ndgrid(xplane, yplane, zplane); %2D cartesian grid of projection plane in plane voxels
  siz        = size(squeeze(X2));
  pos        = [X2(:) Y2(:) Z2(:)]; clear X2 Y2 Z2;
  
  % get the transformation matrix from plotting space to voxel space
  T1 = inv(transform);
  
  % get the transformation matrix to get the projection plane at the right location and orientation into plotting space.
  T2  = [x(:) y(:) ori(:) loc(:); 0 0 0 1];
  
  if 0
    % this is for debugging
    ft_plot_mesh(warp_apply(T2, pos))
    ft_plot_mesh(corner_head)
    axis on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
  end
  
  % get the transformation matrix from projection plane to voxel space
  M = T1*T2;
  
  % get the positions of the pixels of the desires plane in voxel space
  pos = warp_apply(M, pos);
  
  Xi         = reshape(pos(:,1), siz);
  Yi         = reshape(pos(:,2), siz);
  Zi         = reshape(pos(:,3), siz);
  V          = interpn(X,Y,Z,dat,Xi,Yi,Zi,interpmethod);
  [V,Xi,Yi,Zi] = tight(V,Xi,Yi,Zi);
  siz        = size(Xi);
  
  if domask,
    Vmask    = tight(interpn(X,Y,Z,mask,Xi,Yi,Zi,interpmethod));
  end
  
else
  %-------cut a slice without interpolation
  [x, y] = projplane(ori);
  T2     = [x(:) y(:) ori(:) loc(:); 0 0 0 1];
  
  if all(ori==[1 0 0]), xplane = loc(1);   yplane = 1:dim(2); zplane = 1:dim(3); end
  if all(ori==[0 1 0]), xplane = 1:dim(1); yplane = loc(2); zplane = 1:dim(3);   end
  if all(ori==[0 0 1]), xplane = 1:dim(1); yplane = 1:dim(2); zplane = loc(3);   end
  
  [Xi,Yi,Zi] = ndgrid(xplane, yplane, zplane);
  siz        = size(squeeze(Xi));
  Xi         = reshape(Xi, siz);
  Yi         = reshape(Yi, siz);
  Zi         = reshape(Zi, siz);
  V          = reshape(dat(xplane, yplane, zplane), siz);
  if domask,
    Vmask    = reshape(mask(xplane, yplane, zplane), siz);
  end
  
end

% get positions of the plane in plotting space
posh = warp_apply(transform, [Xi(:) Yi(:) Zi(:)], 'homogeneous', 1e-8);
Xh   = reshape(posh(:,1), siz);
Yh   = reshape(posh(:,2), siz);
Zh   = reshape(posh(:,3), siz);

if isempty(cmap),
  %treat as gray value: scale and convert to rgb
  if doscale
    dmin = min(dat(:));
    dmax = max(dat(:));
    V    = (V-dmin)./(dmax-dmin);
  end
  V(isnan(V)) = 0;
  clear dmin dmax;
  % convert anatomy into RGB values
  V = cat(3, V, V, V);
end

if isempty(h),
  h = surface(Xh, Yh, Zh, V);
else
  set(h, 'Cdata', V);
end

set(h, 'linestyle', 'none');
if domask,
  set(h, 'FaceAlpha', 'flat');
  set(h, 'AlphaDataMapping', 'scaled');
  set(h, 'AlphaData', Vmask);
  if ~isempty(opacitylim)
    alim(opacitylim)
  end
end

if ~isempty(cmap)
  colormap(cmap);
end

if ~isempty(clim)
  caxis(clim);
end

% store for future reference
previous_dim  = dim;
previous_mask = mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y] = projplane(z)
[u,s,v] = svd([eye(3) z(:)]);
x = u(:,2)';
y = u(:,3)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,Xi,Yi,Zi] = tight(V,Xi,Yi,Zi)
% cut off the nans at the edges
x = sum(~isfinite(V),1)<size(V,1);
y = sum(~isfinite(V),2)<size(V,2);
V = V(y,x);
if nargin>1, Xi = Xi(y,x); end
if nargin>2, Yi = Yi(y,x); end
if nargin>3, Zi = Zi(y,x); end
