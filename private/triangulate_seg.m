function [pnt, tri] = triangulate_seg(seg, npnt, origin)

% TRIANGULATE_SEG constructs a triangulation of the outer surface of a
% segmented volume. It starts at the center of the volume and projects the
% vertices of an evenly triangulated sphere onto the outer surface. The
% resulting surface is star-shaped from the origin of the sphere.
%
% Use as
%   [pnt, tri] = triangulate_seg(seg, npnt, origin)
%
% Input arguments:
%  seg    = 3D-matrix (boolean) containing segmented volume. If not boolean
%           seg = seg(~=0);
%  npnt   = requested number of vertices
%  origin = 1x3 vector specifying the location of the origin of the sphere
%           in voxel indices. This argument is optional. If undefined, the
%           origin of the sphere will be in the centre of the volume.
%
% Output arguments:
%  pnt = Nx3 matrix of vertex locations
%  tri = Mx3 matrix of triangles
%
% Seg will be checked for holes, and filled if necessary. Also, seg will be
% checked to consist of a single boolean blob. If not, only the outer surface
% of the largest will be triangulated. SPM is used for both the filling and
% checking for multiple blobs.
%
% See also KSPHERE

% Copyright (C) 2005-2012, Robert Oostenveld
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

% impose it to be boolean
seg = (seg~=0);
dim = size(seg);
len = ceil(sqrt(sum(dim.^2))/2);

if ~any(seg(:))
  ft_error('the segmentation is empty')
end

% define the origin if it is not provided in the input arguments
if nargin<3
  origin(1) = dim(1)/2;
  origin(2) = dim(2)/2;
  origin(3) = dim(3)/2;
end

% ensure that the seg consists of only one filled blob.
% if not filled: throw a warning and fill
% if more than one blob: throw a warning and use the biggest

% look for holes
seg = volumefillholes(seg);

% ensure that SPM is available, needed for spm_bwlabel
ft_hastoolbox('spm8up', 3) || ft_hastoolbox('spm2', 1);

% look for >1 blob
[lab, num] = spm_bwlabel(double(seg), 26);
if num>1,
  ft_warning('the segmented volume consists of more than one compartment, using only the biggest one for the segmentation');

  for k = 1:num
    n(k) = sum(lab(:)==k);
  end
  [m,ix] = max(n);
  seg(lab~=ix) = false;
end

% start with a unit sphere with evenly distributed vertices
[pnt, tri] = mesh_sphere(npnt, 'ksphere');

ishollow = false;

for i=1:npnt
  % construct a sampled line from the center of the volume outward into the direction of the vertex
  lin = (0:0.5:len)' * pnt(i,:);
  lin(:,1) = lin(:,1) + origin(1);
  lin(:,2) = lin(:,2) + origin(2);
  lin(:,3) = lin(:,3) + origin(3);
  % round the sampled line towards the nearest voxel indices, which allows
  % a quick nearest-neighbour interpolation/lookup
  lin = round(lin);
  % exclude indices that do not lie within the volume
  sel = lin(:,1)<1 | lin(:,1)>dim(1) | ...
        lin(:,2)<1 | lin(:,2)>dim(2) | ...
        lin(:,3)<1 | lin(:,3)>dim(3);
  lin = lin(~sel,:);
  sel = sub2ind(dim, lin(:,1), lin(:,2), lin(:,3));

  % interpolate the segmented volume along the sampled line
  int = seg(sel);

  % the value along the line is expected to be 1 at first and then drop to 0
  % anything else suggests that the segmentation is hollow
  ishollow = any(diff(int)==1);

  % find the last sample along the line that is inside the segmentation
  sel = find(int, 1, 'last');
  % this is a problem if sel is empty. If so, use the edge of the volume
  if ~isempty(sel)
    % take the last point inside and average with the first point outside
    pnt(i,:) = lin(sel,:);
  else
    % take the edge
    pnt(i,:) = lin(end,:);
  end
end

if ishollow
  % this should not have hapened, especially not after filling the holes
  ft_warning('the segmentation is not star-shaped, please check the surface mesh');
end

% undo the shift of the origin from where the projection is done
% pnt(:,1) = pnt(:,1) - origin(1);
% pnt(:,2) = pnt(:,2) - origin(2);
% pnt(:,3) = pnt(:,3) - origin(3);

% fast unconditional re-implementation of the standard MATLAB function
function [s] = sub2ind(dim, i, j, k)
s = i + (j-1)*dim(1) + (k-1)*dim(1)*dim(2);
