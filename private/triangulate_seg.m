function [pnt, tri] = triangulate_seg(seg, npnt, ori);

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
% $Log: triangulate_seg.m,v $
% Revision 1.4  2009/01/19 12:20:04  roboos
% incorporated suggestion by Jon Iversen to fix bug in case sel=[]
%
% Revision 1.3  2006/07/26 11:05:58  roboos
% use find('last') for matlab 7 and higher, and regular find for older matlab versions
%
% Revision 1.2  2006/04/03 10:39:24  roboos
% added origin of the projection towards the surface as third (optional) input argument
%
% Revision 1.1  2005/11/03 11:12:32  roboos
% new implementation, using a projection of a ksphere triangulation from the center of the volume
%

seg = (seg~=0);
dim = size(seg);
len = ceil(sqrt(sum(dim.^2))/2);

if nargin<3
  ori(1) = dim(1)/2;
  ori(2) = dim(2)/2;
  ori(3) = dim(3)/2;
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

