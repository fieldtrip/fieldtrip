function [pnt, tri] = cylinder(Naz, Nel)

% CYLINDER creates a triangulated cylinder
% 
% Use as
%   [pnt, tri] = cylinder(Naz, Nel)

% Copyright (C) 2002, Robert Oostenveld
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

az = (2*pi*(0:(Naz-1))/Naz)';
el = (linspace(-1,1,Nel))';

pnt = zeros(Naz*Nel+2, 3);
tri = zeros(2*Naz*(Nel-1)+2*Naz,3);

for i=1:Nel
  pnt(((i-1)*Naz+1):(i*Naz), :) = [cos(az) sin(az) el(i)*ones(Naz,1)];
end

pnt(Naz*Nel+1,:) = [0 0 -1];
pnt(Naz*Nel+2,:) = [0 0  1];

ntri = 1;
for i=1:(Nel-1)
  for j=1:(Naz-1)
    tri(ntri, :) = [(i-1)*Naz+j (i-1)*Naz+j+1 i*Naz+j+1];
    ntri = ntri+1;
    tri(ntri, :) = [(i-1)*Naz+j i*Naz+j+1 i*Naz+j];
    ntri = ntri+1;
  end
  tri(ntri, :) = [i*Naz (i-1)*Naz+1 i*Naz+1];
  ntri = ntri+1;
  tri(ntri, :) = [i*Naz i*Naz+1 (i+1)*Naz];
  ntri = ntri+1;
end

% connect the bottom point to the cylinder edge
for i=1:(Naz-1)
  tri(ntri, :) = [Naz*Nel+1 i+1 i];
  ntri = ntri+1;
end
tri(ntri, :) = [Naz*Nel+1 1 Naz];
ntri = ntri+1;

% connect the top point to the cylinder edge
for i=1:(Naz-1)
  tri(ntri, :) = [Naz*Nel+2 Naz*(Nel-1)+i Naz*(Nel-1)+i+1];
  ntri = ntri+1;
end
tri(ntri, :) = [Naz*Nel+2 Naz*Nel Naz*(Nel-1)+1];

