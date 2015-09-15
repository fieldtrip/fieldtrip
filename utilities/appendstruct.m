function s1 = appendstruct(s1, s2);

% APPENDSTRUCT appends a structure to a struct-array. It also works if the initial
% structure is an empty structure or an empty double array.
%
% Use as
%   a = appendstruct(a, b)
% which 
%
% See also PRINTSTRUCT, COPYFIELDS, KEEPFIELDS, REMOVEFIELDS

% Copyright (C) 2015, Robert Oostenveld
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

assert(isstruct(s1) || isempty(s1), 'input argument 1 should be empty or a structure');
assert(isstruct(s2), 'input argument 2 should be a structure');

if isempty(s1)
  s1 = s2;
elseif isstruc
  s1(end+1) = s2;
end
