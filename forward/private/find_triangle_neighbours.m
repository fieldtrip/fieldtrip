function [nb] = find_triangle_neighbours(pnt, tri)

% FIND_TRIANGLE_NEIGHBOURS determines the three neighbours for each triangle
% in a mesh. It returns NaN's if the triangle does not have a neighbour on 
% that particular side.
% 
% [nb] = find_triangle_neighbours(pnt, tri)

% Copyright (C) 2003, Robert Oostenveld
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

npnt = size(pnt,1);
ntri = size(tri,1);

% each triangle has maximally three neighbours, assuming that the 
% surface mesh is not degenerate
nb = nan(size(tri));

% for i=1:ntri
%   for j=setdiff(1:ntri, i)
%     if length(intersect(tri(i,[1 2]), tri(j,:)))==2
%       nb(i,1) = j;
%       continue;
%     end
%     if length(intersect(tri(i,[2 3]), tri(j,:)))==2
%       nb(i,2) = j;
%       continue;
%     end
%     if length(intersect(tri(i,[3 1]), tri(j,:)))==2
%       nb(i,3) = j;
%       continue;
%     end
%   end
%   if all(~isnan(nb(i,:)))
%     continue;
%   end
% end

for i=1:ntri
  % find all neighbouring triangles
  tmp1 = (tri==tri(i,1));
  tmp2 = (tri==tri(i,2));
  tmp3 = (tri==tri(i,3));
  tmp  = (tmp1|tmp2|tmp3);
  sel = find(sum(tmp,2)==2);

  % ensure that each neighbour is assigned to the proper edge
  if length(sel)>3
    error('more than three neighbours found for triangle %d', i);
  else
    for j=1:length(sel)
      if isempty(setdiff(intersect(tri(i,:), tri(sel(j),:)), tri(i,[1 2])))
        nb(i,1) = sel(j);
      elseif isempty(setdiff(intersect(tri(i,:), tri(sel(j),:)), tri(i,[2 3])))
        nb(i,2) = sel(j);
      elseif isempty(setdiff(intersect(tri(i,:), tri(sel(j),:)), tri(i,[3 1])))
        nb(i,3) = sel(j);
      end
    end
  end
end
