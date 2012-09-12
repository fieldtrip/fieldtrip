function [pnt, tri] = triangulate_seg(seg, npnt, ori)

% TRIANGULATE_SEG constructs a triangulation of the outer surface of a
% segmented volume. It starts at the center of the volume and projects the
% vertices of an evenly triangulated sphere onto the outer surface. The
% resulting surface is star-shaped from the center of the volume.
%
% Use as
%   [pnt, tri] = triangulate_seg(seg, npnt)
%
% See also KSPHERE

% Copyright (C) 2005, Robert Oostenveld
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

seg = (seg~=0);
dim = size(seg);
len = ceil(sqrt(sum(dim.^2))/2);

if nargin<3
  ori(1) = dim(1)/2;
  ori(2) = dim(2)/2;
  ori(3) = dim(3)/2;
end

% ensure that the seg consists of only one blob, if not throw a warning and
% use the biggest
ft_hastoolbox('SPM8', 1);
[lab, num] = spm_bwlabel(double(seg), 26);
if num>1,
  warning('the segmented volume consists of more than one compartment, using only the biggest one for the segmentation');
  
  for k = 1:num
    n(k) = sum(lab(:)==k);
  end
  [m,ix] = max(n);
  seg(lab~=ix) = false;
end

% start with a unit sphere with evenly distributed vertices
[pnt, tri] = ksphere(npnt);

for i=1:npnt
  % construct a sampled line from the center of the volume outward into the direction of the vertex
  lin = (0:0.5:len)' * pnt(i,:);
  lin(:,1) = lin(:,1) + ori(1);
  lin(:,2) = lin(:,2) + ori(2);
  lin(:,3) = lin(:,3) + ori(3);
  % round the sampled line towards the nearest voxel indices, which allows
  % a quick nearest-neighbour interpolation/lookup
  lin = round(lin);
  % exclude indices that do not ly within the volume
  sel = lin(:,1)<1 | lin(:,1)>dim(1) | ...
        lin(:,2)<1 | lin(:,2)>dim(2) | ...
        lin(:,3)<1 | lin(:,3)>dim(3);
  lin = lin(~sel,:);
  sel = sub2ind(dim, lin(:,1), lin(:,2), lin(:,3));
  % interpolate the segmented volume along the sampled line
  int = seg(sel);
  % find the last sample along the line that is part of the segmentation
  try
    % for matlab 7 and higher
    sel = find(int, 1, 'last');
  catch
    % for older matlab versions
    sel = find(int);
    sel = sel(end);
  end
  % this is a problem if sel is empty.  If so, use the edge of the volume
  if ~isempty(sel),
      pnt(i,:) = lin(sel,:);
  else
      pnt(i,:) = lin(end,:);
  end
end

% undo the shift of the origin from where the projection is done
% pnt(:,1) = pnt(:,1) - ori(1);
% pnt(:,2) = pnt(:,2) - ori(2);
% pnt(:,3) = pnt(:,3) - ori(3);

% fast unconditional re-implementation of the standard Matlab function
function [s] = sub2ind(dim, i, j, k)
s = i + (j-1)*dim(1) + (k-1)*dim(1)*dim(2);

