function [y] = randsample(x, k)

% RANDSAMPLE Random sample, with or without replacement. This is a drop-in
% replacement for the matlab funnction with the same name. Not all options
% are supported, an error will be issued if needed.
%
% Y = RANDSAMPLE(N,K) returns Y as a 1-by-K vector of values sampled
% uniformly at random, without replacement, from the integers 1:N.
%
% Y = RANDSAMPLE(POPULATION,K) returns K values sampled uniformly at
% random, without replacement, from the values in the vector POPULATION.
%
% See also RAND, RANDPERM.

% Copyright (C) 2007, Robert Oostenveld
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

if nargin>2
  error('only two input variables are supported');
end

if length(x)==1 && isnumeric(x)
  x = 1:x;
end

sel = ceil(rand(1,k)*length(x));
y = x(sel);
