function [vol] = ama2vol(ama)

% AMA2VOL
%
% Use as
%   vol = ama2vol(ama)

% Copyright (C) 2008, Robert Oostenveld
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

vol  = [];
ngeo = length(ama.geo);
for i=1:ngeo
  vol.bnd(i).pos = ama.geo(i).pos;
  vol.bnd(i).tri = ama.geo(i).tri;
  vol.cond(i) = ama.geo(i).sigmam;
end
vol.mat = ama.bi;
npos = size(vol.mat,2);
if size(vol.mat,1)<npos
  vol.mat(npos, npos) = 0;    % it should be a square matrix
end
vol.mat  = vol.mat;
vol.type = 'dipoli';

