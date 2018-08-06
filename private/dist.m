function [d] = dist(x)

% DIST computes the euclidian distance between the columns of the input matrix
%
% Use as
%   [d] = dist(x)
% where x is for example Nx3 for points in 3D space.
%
% This function serves as a replacement for the dist function in the Neural
% Networks toolbox.

% Copyright (C) 2005, Robert Oostenveld
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

n = size(x,2);
d = zeros(n,n);
for i=1:n
  for j=(i+1):n
    d(i,j) = sqrt(sum((x(:,i)-x(:,j)).^2));
    d(j,i) = d(i,j);
  end
end

