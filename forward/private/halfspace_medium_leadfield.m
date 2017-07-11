function [lf] = halfspace_medium_leadfield(rd, elc, vol)

% HALFSPACE_MEDIUM_LEADFIELD calculate the halfspace medium leadfield
% on positions pnt for a dipole at position rd and conductivity cond
% The halfspace solution requires a plane dividing a conductive zone of
% conductivity cond, from a non coductive zone (cond = 0)
%       
% [lf] = halfspace_medium_leadfield(rd, elc, cond)

% Copyright (C) 2011, Cristiano Micheli and Robert Oostenveld
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
  Ndipoles = prod(siz)/3;
  rd = rd(:)'; % ensure that it is a row vector
elseif siz(2)==3
  % positions are specified as a Nx3 matrix -> reformat to a single vector
  Ndipoles = siz(1);
  rd = rd';
  rd = rd(:)'; % ensure that it is a row vector
else
  ft_error('incorrect specification of dipole locations');
end

Nelc     = size(elc,1);
lf       = zeros(Nelc,3*Ndipoles);

for i=1:Ndipoles
  % this is the position of dipole "i"
  dip1 = rd((1:3) + 3*(i-1));
  
  % distances electrodes - dipole
  r1 = elc - ones(Nelc,1) * dip1;
  
  % Method of mirror dipoles:
  % Defines the position of mirror dipoles being symmetric to the plane
  dip2 = get_mirror_pos(dip,vol);
  
  % distances electrodes - mirror dipole
  r2 = elc - ones(Nelc,1) * dip2;
  
  % denominator
  R1 =  (4*pi*vol.cond) * (sum(r1' .^2 ) .^ 1.5)';
  % denominator, mirror term
  R2 = -(4*pi*vol.cond) * (sum(r2' .^2 ) .^ 1.5)';
  
  % condition of dipoles falling in the non conductive halfspace  
  condition = get_dip_halfspace(dip1,vol);
  
  invacuum = acos(dot(ori,(P-pnt)./norm(P-pnt))) < pi/2;
  
  if invacuum
    ft_warning('dipole lies on the vacuum side of the plane');
    lf(:,(1:3) + 3*(i-1)) = NaN(Nelc,3);
  elseif any(R1)==0
    ft_warning('dipole coincides with one of the electrodes');
    lf(:,(1:3) + 3*(i-1)) = NaN(Nelc,3);
  else
    lf(:,(1:3) + 3*(i-1)) = (r ./ [R1 R1 R1]) + (rp ./ [R2 R2 R2]);
  end
end
