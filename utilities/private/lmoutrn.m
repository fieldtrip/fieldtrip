function [la, mu, dist, proj] = lmoutrn(v1, v2, v3, r)

% LMOUTRN computes the la/mu parameters of a point projected to triangles
%
% Use as
%   [la, mu, dist, proj] = lmoutrn(v1, v2, v3, r)
% where v1, v2 and v3 are Nx3 matrices with vertex positions of the triangles, 
% and r is the point that is projected onto the planes spanned by the vertices
% This is a vectorized version of Robert's lmoutrn function and is
% generally faster than a for-loop around the mex-file. It also returns the 
% projection of the point r onto the planes of the triangles, and the signed
% distance to the triangles. The sign of the distance is negative if the point
% lies closer to the average across all vertices and the triangle under consideration.

% Copyright (C) 2012, Jan-Mathijs Schoffelen
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

if size(r,1)==1 && size(v1,1)>1
  r = repmat(r, [size(v1,1), 1]);
end

% compute la/mu parameters
vec0 = r  - v1;
vec1 = v2 - v1;
%vec2 = v3 - v2;
vec3 = v3 - v1;
origin = repmat(mean([v1;v2;v3]), [size(v1,1), 1]);

tmp(:,1,:) = vec1';
tmp(:,2,:) = vec3';
tmp   = pinvNx2(tmp);
la    = sum(vec0'.*shiftdim(tmp(1,:,:))).'; % tmp*vec0';
mu    = sum(vec0'.*shiftdim(tmp(2,:,:))).';

% determine the projection onto the plane of the triangle
proj  = v1 + [la la la].*vec1 + [mu mu mu].*vec3;

% determine the signed distance from the original point to its projection
% where the sign is negative if the original point is closer to the origin 
origin_r    = sum((r    - origin).^2,2);
origin_proj = sum((proj - origin).^2,2);

dist = sqrt(sum((r-proj).^2,2)).*sign(origin_r-origin_proj);

