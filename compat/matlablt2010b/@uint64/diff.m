function y = diff(x)

% DIFF Difference and approximate derivative.
%
% DIFF(X), for a vector X, is [X(2)-X(1)  X(3)-X(2) ... X(n)-X(n-1)].
% DIFF(X), for a matrix X, is the matrix of row differences, [X(2:n,:) - X(1:n-1,:)].
% DIFF(X), for an N-D array X, is the difference along the first non-singleton dimension of X.

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if nargin>1
  error('this implementation is only supported with one input argument');
end

siz = size(x);
if numel(siz)>2
  error('this implementation is only supported with vector or matrix input');
end

if siz(1)==1
  % derivative along the second dimension
  y = x(:,1:end-1);
  for i=1:(siz(2)-1)
    y(:,i) = x(:,end) - y(:,i);
  end

elseif siz(2)==1
  % derivative along the first dimension
  y = x(1:end-1,:);
  for i=1:(siz(1)-1)
    y(i,:) = x(end,:) - y(i,:);
  end

else
  % derivative along the first dimension
  y = x(1:end-1,:);
  for i=1:(siz(1)-1)
    y(i,:) = x(end,:) - y(i,:);
  end
end

