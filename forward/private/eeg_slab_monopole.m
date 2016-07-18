function [lf] = eeg_slab_monopole(rd, elc, vol)

% EEG_SLAB_MONOPOLE calculate the strip medium leadfield
% on positions pnt for a monopole at position rd and conductivity cond
% The halfspace solution requires a plane dividing a conductive zone of
% conductivity cond, from a non coductive zone (cond = 0)
%       
% [lf] = eeg_slab_monopole(rd, elc, cond)
%
% Implemented from Malmivuo J, Plonsey R, Bioelectromagnetism (1993)
% http://www.bem.fi/book/index.htm

% Copyright (C) 2011, Cristiano Micheli
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

siz = size(rd);
if any(siz==1)
  % positions are specified as a single vector
  Npoles = prod(siz)/3;
  rd = rd(:)'; % ensure that it is a row vector
elseif siz(2)==3
  % positions are specified as a Nx3 matrix -> reformat to a single vector
  Npoles = siz(1);
  rd = rd';
  rd = rd(:)'; % ensure that it is a row vector
else
  error('incorrect specification of pole locations');
end

Nelc     = size(elc,1);
lf       = zeros(Nelc,Npoles); 

for i=1:Npoles
  % this is the position of dipole "i"
  pole1 = rd((1:3) + 3*(i-1));
  
  % distances electrodes - corrent poles
  r1 = elc - ones(Nelc,1) * pole1;
  
  % Method of mirror charges:
  % Defines the position of mirror charge being symmetric to the plane
  [pole2,pole3,pole4] = get_mirror_pos(pole1,vol);
  
  % distances electrodes - mirror charge
  r2 = elc - ones(Nelc,1) * pole2;
  r3 = elc - ones(Nelc,1) * pole3;
  r4 = elc - ones(Nelc,1) * pole4;
  
  % denominator
  R1 =  (4*pi*vol.cond) * sqrt(sum(r1' .^2 ) )';
  % denominator, mirror term
  R2 = -(4*pi*vol.cond) * sqrt(sum(r2' .^2 ) )';
  % denominator, mirror term of P1, plane 2
  R3 = -(4*pi*vol.cond) * sqrt(sum(r3' .^2 ) )';  
  % denominator, mirror term of P2, plane 2
  R4 =  (4*pi*vol.cond) * sqrt(sum(r4' .^2 ) )';    
  
  % condition of poles falling in the non conductive halfspace    
  instrip1 = acos(dot(vol.ori1,(pole1-vol.pnt1)./norm(pole1-vol.pnt1))) > pi/2;
  instrip2 = acos(dot(vol.ori2,(pole1-vol.pnt2)./norm(pole1-vol.pnt2))) > pi/2;
  invacuum = ~(instrip1&instrip2);
  
  if invacuum
    warning('a pole lies on the vacuum side of the plane');
    lf(:,i) = NaN(Nelc,1);
  elseif any(R1)==0 
    warning('a pole coincides with one of the electrodes');
    lf(:,i) = NaN(Nelc,1);
  else
    lf(:,i) = (1 ./ R1) + (1 ./ R2) + (1 ./ R3);% + (1 ./ R4);
  end
end

function [P2,P3,P4] = get_mirror_pos(P1,vol)
% calculates the position of a point symmetric to the pole, with respect to plane1
% and two points symmetric to the last ones, with respect to plane2

P2 = []; P3 = []; P4 = [];

% define the planes
pnt1 = vol.pnt1;
ori1 = vol.ori1;
pnt2 = vol.pnt2;
ori2 = vol.ori2;

if abs(dot(P1-pnt1,ori1))<eps || abs(dot(P1-pnt2,ori2))<eps
  warning(sprintf ('point %f %f %f lies on the plane',P1(1),P1(2),P1(3)))
  P2 = P1;
else
  
  % define the planes
  plane1 = def_plane(pnt1,ori1);
  plane2 = def_plane(pnt2,ori2);
  
  % distance plane1-point P1
  d = abs(dot(ori1, plane1(:,1:3)-P1(:,1:3), 2));
  % symmetric point
  P2 = P1 + 2*d*ori1;
  
  % distance plane2-point P1
  d = abs(dot(ori2, plane2(:,1:3)-P1(:,1:3), 2));  
  % symmetric point
  P3 = P1 + 2*d*ori2;  
  
  % distance plane2-point P2
  d = abs(dot(ori2, plane2(:,1:3)-P2(:,1:3), 2));  
  % symmetric point
  P4 = P2 + 2*d*ori2;    
end

function plane = def_plane(pnt,ori)
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
