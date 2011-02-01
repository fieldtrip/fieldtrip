function [lf] = halfspace_medium_leadfield(rd, pnt, vol);

% HALFSPACE_MEDIUM_LEADFIELD calculate the halfspace medium leadfield
% on positions pnt for a dipole at position rd and conductivity cond
% The halfspace solution requires a plane dividing a conductive zone of
% conductivity cond, from a non coductive zone (cond = 0)
%       
% [lf] = halfspace_medium_leadfield(rd, pnt, cond)

% Copyright (C) 2011, Cristiano Micheli and Robert Oostenveld
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
% $Id: halfspace_medium_leadfield.m $

siz = size(rd);
if any(siz==1)
  % positions are specified as a single vector
  Ndipoles = prod(siz)/3;
elseif siz(2)==3
  % positions are specified as a Nx3 matrix -> reformat to a single vector
  Ndipoles = siz(1);
  rd = rd';
  rd = rd(:);
else
  error('incorrect specification of dipole locations');
end

Npnt     = size(pnt,1);
lf       = zeros(Npnt,3*Ndipoles);
s1       = size(rd);

if s1(1)>s1(2)
  % make sure that the dipole position is a row-vector
  rd = rd';
end

for i=1:Ndipoles
  % distances sensors - dipole
  r = pnt - ones(Npnt,1) * rd((1:3) + 3*(i-1));
  
  % Method of mirror dipoles:
  % Defines the position of mirror dipoles being symmetric to the plane
  rdp = get_mirror_pos(rd((1:3) + 3*(i-1)),vol);
  
  % distances sensors - mirror dipole
  rp = pnt - ones(Npnt,1) * rdp;
  
  % denominator
  R1 =  (4*pi*vol.cond) * (sum(r' .^2 ) .^ 1.5)';
  % denominator, mirror term
  R2 = -(4*pi*vol.cond) * (sum(rp' .^2 ) .^ 1.5)';
  
  % condition of dipoles falling in the non conductive halfspace  
  condition = get_dip_halfspace(rd((1:3) + 3*(i-1)),vol);
  
  if any(R1)==0
    warning('dipole lies on boundary of volume model');
  elseif condition
    lf(:,(1:3) + 3*(i-1)) = zeros(size(r,1),3);
  else
    lf(:,(1:3) + 3*(i-1)) = (r ./ [R1 R1 R1]) + (rp ./ [R2 R2 R2]);
  end
end
