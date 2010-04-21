function [x,mx,sx] = standardise(x,dim)

% X = NANSTANDARDISE(X, DIM) computes the zscore of a matrix along dimension 
% dim, taking nans into account
% has similar functionality as the stats-toolbox's zscore function

% Copyright (C) 2010, Jan-Mathijs Schoffelen
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

if nargin == 1, dim = find(size(x)>1,1,'first'); end

siz    = size(x);
n      = sum(~isnan(x),dim);
x(isnan(x)) = 0;
repsiz = siz;
repsiz(setdiff(1:numel(siz), dim)) = 1;
ressiz      = [siz 1];
if dim>1,
  ressiz(dim) = [];
else
  ressiz(dim) = 1;
end
mx     = sum(x,dim)./n;
x      = x - repmat(mx, repsiz);
%mx     = reshape(mx, ressiz);
sx     = sqrt(sum(x.^2,dim)./n);
x      = x ./repmat(sx, repsiz);
%sx     = reshape(sx, ressiz);
