function [lf] = inf_medium_leadfield(rd, pnt, cond)

% INF_MEDIUM_LEADFIELD calculate the infinite medium leadfield
% on positions pnt for dipole position R and conductivity cond
%       
% [lf] = inf_medium_leadfield(R, pnt, cond)

% Copyright (C) 1998, Robert Oostenveld
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
elseif siz(2)==3
  % positions are specified as a Nx3 matrix -> reformat to a single vector
  Ndipoles = siz(1);
  rd = rd';
  rd = rd(:);
else
  ft_error('incorrect specification of dipole locations');
end

Npnt     = size(pnt,1);
lf       = zeros(Npnt,3*Ndipoles);
s1       = size(rd);

if s1(1)>s1(2)
  % make sure that the dipole position is a row-vector
  rd = rd';
end

for i=1:Ndipoles
  r = pnt - ones(Npnt,1) * rd((1:3) + 3*(i-1));
  R = (4*pi*cond) * (sum(r' .^2 ) .^ 1.5)';
  if any(R==0)
    ft_warning('dipole lies on boundary of volume model');
  end
  lf(:,(1:3) + 3*(i-1)) = r ./ [R R R];
end

