function [X, Y, Z, pos1, tri1, pos2, tri2] = intersect_plane(pos, tri, v1, v2, v3)

% INTERSECT_PLANE intersection between a triangulated surface mesh and a plane. It
% returns the coordinates of the begin- and endpoints of the line segments that
% together form the contour of the intersection.
%
% Use as
%   [X, Y, Z] = intersect_plane(pos, tri, v1, v2, v3)
%
% where the intersecting plane is spanned by the vertices v1, v2, v3 and the return
% values are the X, Y and Z coordinates of the begin- and endpoints for all line
% segments.

% Copyright (C) 2002-2024, Robert Oostenveld
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

if ~isa(pos, 'double'), pos = double(pos); end % low level mex files require double precision input

% determine on which side of the plane all vertices are
side = ptriside(v1, v2, v3, pos);

% remember the original mesh description
original = tri;

if size(tri,2)==4
  % the input mesh consists of tetraheders

  % only continue with the tetraheders that have vertices on both sides of the plane
  indx = abs(sum(side(tri),2))~=4;
  tri = tri(indx,:);

  if ~isempty(tri)
    % construct and concatenate the four triangular edges of each tetraheder
    part1 = tri(:,[1 2 3]);
    part2 = tri(:,[1 2 4]);
    part3 = tri(:,[2 3 4]);
    part4 = tri(:,[3 1 4]);
    tri = cat(1, part1, part2, part3, part4);
  end
elseif size(tri,2)==8
  % the input mesh consists of hexaheders

  % only continue with the hexaheders that have vertices on both sides of the plane
  indx = abs(sum(side(tri),2))~=8;
  tri = tri(indx,:);

  if ~isempty(tri)
    cube.hex = [1 2 3 4 5 6 7 8];
    cube.pos = pos(tri(1,:),:);
    cube.tri = convhull(cube.pos(:,1), cube.pos(:,2), cube.pos(:,3));
    % each hexaheder/cube can be described by 12 triangles, 2 for each side
    part1  = tri(:,cube.tri(1,:));
    part2  = tri(:,cube.tri(2,:));
    part3  = tri(:,cube.tri(3,:));
    part4  = tri(:,cube.tri(4,:));
    part5  = tri(:,cube.tri(5,:));
    part6  = tri(:,cube.tri(6,:));
    part7  = tri(:,cube.tri(7,:));
    part8  = tri(:,cube.tri(8,:));
    part9  = tri(:,cube.tri(9,:));
    part10  = tri(:,cube.tri(10,:));
    part11  = tri(:,cube.tri(11,:));
    part12  = tri(:,cube.tri(12,:));
    tri = cat(1, part1, part2, part3, part4, part5, part6, part7, part8, part9, part10, part11, part12);
  end
end

% find the triangles which have vertices on both sides of the plane
indx = find(abs(sum(side(tri),2))~=3);
cnt1 = zeros(length(indx), 3);
cnt2 = zeros(length(indx), 3);

for i=1:length(indx)
  cur = tri(indx(i),:);
  tmp = side(cur);
  l1 = pos(cur(1),:);
  l2 = pos(cur(2),:);
  l3 = pos(cur(3),:);
  if tmp(1)==tmp(2)
    % plane intersects two sides of the triangle
    cnt1(i,:) = ltrisect(v1, v2, v3, l3, l1);
    cnt2(i,:) = ltrisect(v1, v2, v3, l3, l2);
  elseif tmp(1)==tmp(3)
    cnt1(i,:) = ltrisect(v1, v2, v3, l2, l1);
    cnt2(i,:) = ltrisect(v1, v2, v3, l2, l3);
  elseif tmp(2)==tmp(3)
    cnt1(i,:) = ltrisect(v1, v2, v3, l1, l2);
    cnt2(i,:) = ltrisect(v1, v2, v3, l1, l3);
  elseif tmp(1)==0 && tmp(2)==0
    % two vertices of the triangle lie on the plane
    cnt1(i,:) = l1;
    cnt2(i,:) = l2;
  elseif tmp(1)==0 && tmp(3)==0
    cnt1(i,:) = l1;
    cnt2(i,:) = l3;
  elseif tmp(2)==0 && tmp(3)==0
    cnt1(i,:) = l2;
    cnt2(i,:) = l3;
  elseif tmp(1)==0 && tmp(2)~=tmp(3)
    % one vertex of the triangle lies on the plane
    cnt1(i,:) = l1;
    cnt2(i,:) = ltrisect(v1, v2, v3, l2, l3);
  elseif tmp(2)==0 && tmp(3)~=tmp(1)
    cnt1(i,:) = l2;
    cnt2(i,:) = ltrisect(v1, v2, v3, l3, l1);
  elseif tmp(3)==0 && tmp(1)~=tmp(2)
    cnt1(i,:) = l3;
    cnt2(i,:) = ltrisect(v1, v2, v3, l1, l2);
  elseif tmp(1)==0
    cnt1(i,:) = l1;
    cnt2(i,:) = l1;
  elseif tmp(2)==0
    cnt1(i,:) = l2;
    cnt2(i,:) = l2;
  elseif tmp(3)==0
    cnt1(i,:) = l3;
    cnt2(i,:) = l3;
  end
end

X = [cnt1(:,1) cnt2(:,1)];
Y = [cnt1(:,2) cnt2(:,2)];
Z = [cnt1(:,3) cnt2(:,3)];

if nargout>3
  % also output the two meshes on either side of the plane
  [pos1, tri1] = remove_vertices(pos, original, side~= 1);
  [pos2, tri2] = remove_vertices(pos, original, side~=-1);
end
