function [x, mx, sx] = standardise(x, dim, lim)

% STANDARDISE computes the zscore of a matrix along dimension dim
% has similar functionality as the stats-toolbox's zscore function
%
% Use as
%   x = standardise(x, dim)
%
% See also ZSCORE

% Copyright (C) 2009, Jan-Mathijs Schoffelen
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

if nargin == 1,
  dim = find(size(x)>1,1,'first');
end

if nargin == 3,
  error('third input argument is not used');
end

switch dim
case 1
  n  = size(x,dim);
  mx = mean(x,dim);
  x  = x-mx(ones(1,n),:,:,:,:,:,:,:);
  sx = sqrt(sum(abs(x).^2,dim)./n);
  x  = x./sx(ones(1,n),:,:,:,:,:,:,:);
case 2
  n  = size(x,dim);
  mx = mean(x,dim);
  x  = x-mx(:,ones(1,n),:,:,:,:,:,:);
  sx = sqrt(sum(abs(x).^2,dim)./n);
  x  = x./sx(:,ones(1,n),:,:,:,:,:,:);
case 3
  n  = size(x,dim);
  mx = mean(x,dim);
  x  = x-mx(:,:,ones(1,n),:,:,:,:,:);
  sx = sqrt(sum(abs(x).^2,dim)./n);
  x  = x./sx(:,:,ones(1,n),:,:,:,:,:);
case 4
  n  = size(x,dim);
  mx = mean(x,dim);
  x  = x-mx(:,:,:,ones(1,n),:,:,:,:);
  sx = sqrt(sum(abs(x).^2,dim)./n);
  x  = x./sx(:,:,:,ones(1,n),:,:,:,:);
case 5
  n  = size(x,dim);
  mx = mean(x,dim);
  x  = x-mx(:,:,:,:,ones(1,n),:,:,:);
  sx = sqrt(sum(abs(x).^2,dim)./n);
  x  = x./sx(:,:,:,:,ones(1,n),:,:,:);
case 6
  n  = size(x,dim);
  mx = mean(x,dim);
  x  = x-mx(:,:,:,:,:,ones(1,n),:,:);
  sx = sqrt(sum(abs(x).^2,dim)./n);
  x  = x./sx(:,:,:,:,:,ones(1,n),:,:);
case 7
  n  = size(x,dim);
  mx = mean(x,dim);
  x  = x-mx(:,:,:,:,:,:,ones(1,n),:);
  sx = sqrt(sum(abs(x).^2,dim)./n);
  x  = x./sx(:,:,:,:,:,:,ones(1,n),:);
case 8
  n  = size(x,dim);
  mx = mean(x,dim);
  x  = x-mx(:,:,:,:,:,:,:,ones(1,n));
  sx = sqrt(sum(abs(x).^2,dim)./n);
  x  = x./sx(:,:,:,:,:,:,:,ones(1,n));
otherwise
  error('dim too large, standardise currently supports dimensionality up to 8');
end
