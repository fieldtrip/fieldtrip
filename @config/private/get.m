function y = get(x, key, inc)

% GET Return the value of a field in a config object.

% Copyright (C) 2012, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if nargin<3
  % the default is to increment the reference counter, an exception is when
  % a structure has to be recursed to assign a value to a substructure
  inc = true;
end

if isfield(x.value, key)
  y = x.value.(key);
  if inc
    increment(x.reference.(key));
  end
else
  error(sprintf('Reference to non-existent field ''%s''.', key));
end
