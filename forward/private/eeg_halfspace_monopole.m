function [lf] = eeg_halfspace_monopole(rd, elc, vol)

% EEG_HALFSPACE_MONOPOLE calculate the halfspace medium leadfield
% on positions pnt for a monopole at position rd and conductivity cond
% The halfspace solution requires a plane dividing a conductive zone of
% conductivity cond, from a non coductive zone (cond = 0)
%       
% [lf] = eeg_halfspace_monopole(rd, elc, cond)
%
% Implemented from Malmivuo J, Plonsey R, Bioelectromagnetism (1993)
% http://www.bem.fi/book/index.htm

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
% $Id: eeg_halfspace_monopole.m $

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
  error('incorrect specification of dipole locations');
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
  pole2 = get_mirror_pos(pole1,vol);
  
  % distances electrodes - mirror charge
  r2 = elc - ones(Nelc,1) * pole2;
  
  % denominator
  R1 =  (4*pi*vol.cond) * (sum(r1' .^2 ) )';
  % denominator, mirror term
  R2 = -(4*pi*vol.cond) * (sum(r2' .^2 ) )';
  
  % condition of poles falling in the non conductive halfspace    
  invacuum = acos(dot(vol.ori,(pole1-vol.pnt)./norm(pole1-vol.pnt))) < pi/2;
  
  if invacuum
    warning('a pole lies on the vacuum side of the plane');
    lf(:,i) = NaN(Nelc,1);
  elseif any(R1)==0
    warning('a pole coincides with one of the electrodes');
    lf(:,i) = NaN(Nelc,1);
  else
    lf(:,i) = (1 ./ R1) + (1 ./ R2);
  end
end

function P2 = get_mirror_pos(P1,vol)
% calculates the position of a point symmetric to pnt with respect to a plane

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
  d = abs(dot(ori, plane(:,1:3)-P1(:,1:3), 2));

  % symmetric point
  P2 = P1 + 2*d*ori;
end
