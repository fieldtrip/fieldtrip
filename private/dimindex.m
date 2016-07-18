function M = dimindex(A,dim,idx)

% DIMINDEX makes a selection from a multi-dimensional array where the dimension is
% selected by a scalar, not by the place between the brackets.
%
% Use as
%   M = dimindex(A,dim,idx)
%
% The purpose of the function is shown by the following example:
%
% A(:,:,:,23,:,:,...) is the same as dimindex(A,4,23)
% A(2,4,3)            is the same as dimindex(A,[1,2,3],[2,4,3])
% A(4,:,[5:10])       is the same as dimindex(A,[1,3],{4,[5:10]})
%
% See also the function DIMASSIGN

% Copyright (C) 2005, Geerten Kramer
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

if (~iscell(idx))
  if (~any(size(dim)==1) || ~any(size(idx)==1) || ndims(dim)>2 || ndims(idx)>2 || length(dim)~=length(idx))
    error('dim and idx must be both scalars or both vectors of the same size');
  end
  dummi = [];
  for i=1:length(idx)
    dummi{i} = idx(i);
  end
  idx = dummi;
  clear dummi;
end
if (~any(size(dim)==1) || ~any(size(idx)==1) || ndims(dim)>2 || ndims(idx)>2 || length(dim)~=length(idx))
  error('dim and idx must be both scalars or both must have the same length');
end


if (length(dim)>ndims(A)||any(dim>ndims(A)))
  error('dim, or one of its contents are larger than the number of dimentions in A');
end
if (~isequal(unique(dim),sort(dim)))
  error('dim must be unique, every dimention can be addressed only once');
end

N = ndims(A);
for i=1:N
  ref = find(dim==i);
  if (isempty(ref))
    C{i} = ':';
  else
    C{i} = idx{ref};
  end
end
M = A(C{:});

