function is_in_empty = get_dip_halfspace(P,vol)

% GET_DIP_HALFSPACE checks if the dipole is in the empty halfspace and
% returns true if this happens. The normal of the plane points by
% convention to the empty halfspace.
%       
% is_in_empty = get_dip_halfspace(P,vol);

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

is_in_empty = false;
ori = vol.ori;
pnt = vol.pnt;

if acos(dot(ori,(P-pnt)./norm(P-pnt))) < pi/2
  is_in_empty = true;
end
