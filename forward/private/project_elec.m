function [el] = project_elec(elc, pnt, tri)

% PROJECT_ELEC projects electrodes on a triangulated surface
% and returns triangle index, la/mu parameters and distance
%
% [el] = project_elec(elc, pnt, tri)
%
% it returns a Nx4 matrix with [tri, la, mu, dist] for each electrode
%
% See also TRANSFER_ELEC 

% Copyright (C) 1999-2002, Robert Oostenveld
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

Nelc = size(elc,1);
Npnt = size(pnt,1);
Ntri = size(tri,1);
el   = zeros(Nelc, 4);

for i=1:Nelc
  smallest_dist = Inf;

  for j=1:Ntri
    [proj, dist] = ptriproj(pnt(tri(j,1),:), pnt(tri(j,2),:), pnt(tri(j,3),:), elc(i,:), 1);
    if dist<smallest_dist
      % remember the triangle index, distance and la/mu
      [la, mu] = lmoutr(pnt(tri(j,1),:), pnt(tri(j,2),:), pnt(tri(j,3),:), proj);
      smallest_dist = dist; 
      smallest_tri  = j; 
      smallest_la   = la; 
      smallest_mu   = mu; 
    end
  end

  % store the projection for this electrode
  el(i,:) = [smallest_tri smallest_la smallest_mu smallest_dist];
end

