function [pnt3, tri3] = retriangulate(pnt1, tri1, pnt2, tri2, flag)

% RETRIANGULATE projects a triangulation onto another triangulation
% thereby providing a a new triangulation of the old one.
% 
% Use as
%   [pnt, tri] = retriangulate(pnt1, tri1, pnt2, tri2, flag)
% where
%   pnt1, tri1  describe the desired triagulated surface shape
%   pnt2, tri2  describe the triangulation that will be projected
%
% The optional flag can be 0/1/2 and determines whether the center
% of the triangulations should be shifted to the origin before the
% projection is done. The resulting surface will be shifted back to
% its original location.
% 
% flag=0 means no shift (default)
% flag=1 means shifting to the geometrical mean
% flag=2 means shifting to the center of the bounding box
%
% The projection is done from the coordinate system origin (0,0,0). 
%
% See also ICOSAHEDRONxxx, ISOSURFACE, REDUCEPATCH

% Copyright (C) 2003, Robert Oostenveld

if nargin<5
  flag = 0;
end

if flag==0
  center1 = [0 0 0];
  center2 = [0 0 0];
elseif flag==1
  center1 = mean(pnt1,2);
  center2 = mean(pnt2,2);
elseif flag==2
  center1 = (min(pnt1) + max(pnt1))./2;
  center2 = (min(pnt2) + max(pnt2))./2;
end

% shift the center of both surfaces
pnt1(:,1) = pnt1(:,1) - center1(1);
pnt1(:,2) = pnt1(:,2) - center1(2);
pnt1(:,3) = pnt1(:,3) - center1(3);
pnt2(:,1) = pnt2(:,1) - center2(1);
pnt2(:,2) = pnt2(:,2) - center2(2);
pnt2(:,3) = pnt2(:,3) - center2(3);

npnt1 = size(pnt1,1);
npnt2 = size(pnt2,1);

ntri1 = size(tri1,1);
ntri2 = size(tri2,1);

radius1  = sum(pnt1.^2, 2).^0.5;
radius2  = sum(pnt2.^2, 2).^0.5;
cosangle = (pnt1 * pnt2') ./ (radius1 * radius2');
[m, indx] = min(cosangle, [], 2);

pnt3  = zeros(size(pnt2));
count = zeros(size(pnt2,1),1);
for i=1:length(indx)
  count(indx(i))  = count(indx(i)) + 1;
  pnt3(indx(i),:) = pnt3(indx(i),:) + pnt1(i,:);
end

% turn the warning off
ws = warning('off', 'MATLAB:divideByZero');
pnt3 = pnt3 ./ [count count count];
tri3 = tri2;
% revert to previous warning state
warning(ws)

exception=find(count==0);
if ~isempty(exception)
  % warning('treating %d vertices as exception', length(exception));
  [m, indx] = min(cosangle, [], 1);
  for i=exception(:)'
    sel = find_vertex_neighbours(pnt3, tri3, i);
    sel = setdiff(sel, exception);
    if ~isempty(sel)
      pnt3(i,:) = mean(pnt3(sel, :), 1);
    else
      pnt3(i,:) = pnt1(indx(i),:);
    end
  end
end

% shift the center of the new surface back to the original center of surface1
pnt3(:,1) = pnt3(:,1) + center1(1);
pnt3(:,2) = pnt3(:,2) + center1(2);
pnt3(:,3) = pnt3(:,3) + center1(3);

