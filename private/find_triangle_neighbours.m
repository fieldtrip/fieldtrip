function [nb] = find_triangle_neighbours(pnt, dhk)

% FIND_TRIANGLE_NEIGHBOURS determines the three neighbours for each triangle
% in a mesh. It returns NaN's if the triangle does not have a neighbour on 
% that particular side.
% 
% [nb] = find_triangle_neighbours(pnt, dhk)

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: find_triangle_neighbours.m,v $
% Revision 1.2  2003/03/04 21:46:18  roberto
% added CVS log entry and synchronized all copyright labels
%

npnt = size(pnt,1);
ndhk = size(dhk,1);

% each triangle has maximally three neighbours, assuming that the 
% surface mesh is not degenerate
nb = nan * ones(size(dhk));

% for i=1:ndhk
%   for j=setdiff(1:ndhk, i)
%     if length(intersect(dhk(i,[1 2]), dhk(j,:)))==2
%       nb(i,1) = j;
%       continue;
%     end
%     if length(intersect(dhk(i,[2 3]), dhk(j,:)))==2
%       nb(i,2) = j;
%       continue;
%     end
%     if length(intersect(dhk(i,[3 1]), dhk(j,:)))==2
%       nb(i,3) = j;
%       continue;
%     end
%   end
%   if all(~isnan(nb(i,:)))
%     continue;
%   end
% end

for i=1:ndhk
  % find all neighbouring triangles
  tmp1 = (dhk==dhk(i,1));
  tmp2 = (dhk==dhk(i,2));
  tmp3 = (dhk==dhk(i,3));
  tmp  = (tmp1|tmp2|tmp3);
  sel = find(sum(tmp,2)==2);

  % ensure that each neighbour is assigned to the proper edge
  if length(sel)>3
    error(sprintf('more than three neighbours found for triangle %d', i));
  else
    for j=1:length(sel)
      if isempty(setdiff(intersect(dhk(i,:), dhk(sel(j),:)), dhk(i,[1 2])))
        nb(i,1) = sel(j);
      elseif isempty(setdiff(intersect(dhk(i,:), dhk(sel(j),:)), dhk(i,[2 3])))
        nb(i,2) = sel(j);
      elseif isempty(setdiff(intersect(dhk(i,:), dhk(sel(j),:)), dhk(i,[3 1])))
        nb(i,3) = sel(j);
      end
    end
  end
end
