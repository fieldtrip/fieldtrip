function [s] = removefields(s, fields)

% REMOVEFIELDS makes a selection of the fields in a structure
%
% Use as
%   s = removefields(s, fields);
%
% See also KEEPFIELDS, COPYFIELDS

% Copyright (C) 2014, Robert Oostenveld
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

if isempty(s)
   % this prevents problems if s is an empty double, i.e. []
  return
end

if ischar(fields)
  fields = {fields};
elseif ~iscell(fields)
  error('fields input argument must be a cell array of strings or a single string');
end

fields = intersect(fieldnames(s), fields);
for i=1:numel(fields)
  s = rmfield(s, fields{i});
end
