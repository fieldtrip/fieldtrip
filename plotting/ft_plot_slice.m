function [h] = ft_plot_slice(dat, varargin)

% FT_PLOT_SLICE cuts a 2-D slice from a 3-D volume and interpolates
% if necessary
%
% Use as
%   ft_plot_slice(dat, ...)
%
% Additional options should be specified in key-value pairs and can be
%   'transform'    a 4x4 homogeneous transformation matrix specifying the mapping from
%                    voxel space to the coordinate system in which the data are plotted.
%   'location'     a 1x3 vector specifying a point on the plane which will be plotted
%                    the coordinates are expressed in the coordinate system in which the 
%                    data will be plotted. location defines the origin of the plane 
%   'orientation'  a 1x3 vector specifying the direction orthogonal through the plane
%                    which will be plotted.
%
%   'mask'
%   'interplim'
%   'interpmethod'
%   'style'
%   'datmask'

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

% these are for speeding up the plotting on subsequent calls
persistent previous_argin previous_maskimage

transform = keyval('transform',   varargin); if isempty(transform),  transform  = eye(4);    end;
loc       = keyval('location',    varargin); if isempty(loc),        loc        = [0 0 0];   end;
ori       = keyval('orientation', varargin); if isempty(ori),        ori        = [0 0 1];   end;
resolution = keyval('resolution', varargin); if isempty(resolution), resolution = 1;         end
if ~strcmp(class(dat), 'double'),
  origclass = class(dat);
  dat       = cast(dat, 'double');
end

ori = ori./sqrt(sum(ori.^2));
dim = size(dat);

% determine whether interpolation is needed
dointerp = false;
if ~dointerp && sum(sum(transform-eye(4)))~=0, dointerp = true; end
if ~dointerp && ~all(round(loc)==loc),   dointerp = true; end
if ~dointerp && sum(ori)~=1,             dointerp = true; end
if ~dointerp && ~(resolution==round(resolution)), dointerp = true; end

if dointerp
  % get voxel indices
  [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));

  % define 'x' and 'y' axis in projection plane.
  % this is more or less arbitrary
  [x, y] = projplane(ori);
  
  % get transformation matrix from projection space to plotting space.
  T = inv([ [x; y; ori; [0 0 0]] [loc(:); 1]]);
 
  m      = max(dim)./1.5; %1.5 is arbitrary
  xplane = -m:resolution:m;
  yplane = -m:resolution:m;
  zplane = 0;
  [X2,Y2,Z2] = ndgrid(xplane, yplane, zplane); %2D cartesian grid of projection plane
  siz        = size(squeeze(X2));
  pos        = [X2(:) Y2(:) Z2(:)]; clear X2 Y2 Z2;
  pos        = warp_apply(inv(T*transform), pos); % gives the positions of the required plane in voxel space
  Xi         = reshape(pos(:,1), siz);
  Yi         = reshape(pos(:,2), siz);
  Zi         = reshape(pos(:,3), siz);
  V          = interpn(X,Y,Z,dat,Xi,Yi,Zi,'nearest');

else
  if all(ori==[1 0 0]), xplane = loc(1); yplane = 1:dim(2); zplane = 1:dim(3); end
  if all(ori==[0 1 0]), xplane = 1:dim(1); yplane = loc(2); zplane = 1:dim(3); end
  if all(ori==[0 0 1]), xplane = 1:dim(1); yplane = 1:dim(2); zplane = loc(3); end
  [Xi,Yi,Zi] = ndgrid(xplane, yplane, zplane);
  siz        = size(squeeze(Xi));
  Xi         = reshape(Xi, siz);
  Yi         = reshape(Yi, siz);
  Zi         = reshape(Zi, siz);
  V          = dat(xplane, yplane, zplane);

end

% get positions of the plane in head space
posh = warp_apply(transform, [Xi(:) Yi(:) Zi(:)]);
Xh   = reshape(posh(:,1), siz);
Yh   = reshape(posh(:,2), siz);
Zh   = reshape(posh(:,3), siz);

h = surf(Xh, Yh, Zh, V); axis equal
%set(h, 'linestyle', 'none', 'facealpha', 0.9);
set(h, 'linestyle', 'none');

function [x, y] = projplane(z)

[u,s,v] = svd([eye(3) z(:)]);

x = u(:,2)';
y = u(:,3)';
