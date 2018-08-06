function [tra] = transfer_elec(pnt, tri, el); 

% TRANSFER_ELEC is the transfermatrix from vertex to electrode potentials
% using bi-linear interpolation over the triangles
%
% tra = transfer_elec(pnt, tri, el)
%
% the Nx3 matrix el shold contain [tri, la, mu] for each electrode
%
% See also PROJECT_ELEC

% Copyright (C) 1998-2002, Robert Oostenveld
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

Npnt = size(pnt,1);
Ntri = size(tri,1);
Nel  = size(el,1);

tra = zeros(Nel, Npnt);

for i=1:Nel
  tra(i, tri(el(i,1), 1)) = 1 - el(i,2) - el(i,3);
  tra(i, tri(el(i,1), 2)) = el(i,2);
  tra(i, tri(el(i,1), 3)) = el(i,3);
end

