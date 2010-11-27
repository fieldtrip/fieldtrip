function [vol] = ama2vol(ama)

% AMA2VOL
%
% Use as
%   vol = ama2vol(ama)

% Copyright (C) 2008, Robert Oostenveld
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
% $Id: ama2vol.m 945 2010-04-21 17:41:20Z roboos $

vol  = [];
ngeo = length(ama.geo);
for i=1:ngeo
  vol.bnd(i).pnt = ama.geo(i).pnt;
  vol.bnd(i).tri = ama.geo(i).dhk;
  vol.cond(i) = ama.geo(i).sigmam;
end
vol.mat = ama.bi;
npnt = size(vol.mat,2);
if size(vol.mat,1)<npnt
  vol.mat(npnt, npnt) = 0;    % it should be a square matrix
end
vol.mat  = sparse(vol.mat);   % convert to sparse for faster multiplications
vol.type = 'dipoli';

