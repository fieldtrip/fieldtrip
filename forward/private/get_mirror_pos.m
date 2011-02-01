function P2 = get_mirror_pos(P1,vol)

% GET_MIRROR_POS calculates the position of a point symmetric to pnt with respect to a plane
% 
% P1 is a [1x3] point
% vol is the headmodel for halfspace
%
% P2 = get_mirror_pos(P1,vol);

% Copyright (C) 2011, Cristiano Micheli 
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
% $Id: get_mirror_pos.m $

P2 = [];

% define the plane
pnt = vol.pnt;
ori = vol.ori; % already normalized

if abs(dot(P1-pnt,ori))<eps
  warning(sprintf ('point %f %f %f lies in the symmetry plane',P1(1),P1(2),P1(3)))
  P2 = P1;
else
  % define the plane in parametric form
  % define a non colinear vector vc with respect to the plane normal
  vc = [1 0 0];    
  if abs(cross(ori, vc, 2))<eps
      vc = [0 1 0];
  end
  % define plane's direction vectors
  v1 = cross(ori, vc, 2);  v1 = v1/norm(v1);
  v2 = cross(pnt, ori, 2); v2 = v2/norm(v2);
  plane = [pnt v1 v2];
 
  % distance plane-point P1
  d = -dot(ori, plane(:,1:3)-P1(:,1:3), 2);

  % symmetric point
  P2 = P1 + 2*d*ori;
end

