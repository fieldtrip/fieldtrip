function [d] = dist(x, y)

% DIST computes the Euclidian distance between the columns of the input matrix or
% between the rows and columns of two input matrices.
%
% This function serves as a drop-in replacement for the dist function in the Neural
% Networks toolbox.
%
% Use as
%   [d] = dist(x')
% where x is for example an Nx3 matrix with vertices in 3D space, or as
%   [d] = dist(x, y')
% where x and y are Nx3 and Mx3 matrices with vertices in 3D space
%
% See also DSEARCHN, KNNSEARCH

% Copyright (C) 2005-2024, Robert Oostenveld
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

if nargin==1
  n = size(x,2);
  d = zeros(n,n);
  for i=1:n
    for j=(i+1):n
      d(i,j) = sqrt(sum((x(:,i)-x(:,j)).^2));
      d(j,i) = d(i,j);
    end
  end

elseif nargin==2

  n = size(x,1);
  m = size(y,2);
  d = zeros(n,m);

  if m==1
    % do it the efficient way
    x(:,1) = x(:,1) - y(1);
    x(:,2) = x(:,2) - y(2);
    x(:,3) = x(:,3) - y(3);
    d = sqrt(sum(x.^2,2));

  else
    % do it the normal way
    for i=1:n
      for j=1:m
        d(i,j) = sqrt(sum((x(i,:)-y(:,j)').^2));
      end
    end
  end

end % if nargin is 1 or 2
