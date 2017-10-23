function b = renamefields(a, old, new)

% RENAMEFIELDS renames a selection of the fields in a structure
%
% Use as
%   b = renamefields(a, old, new)
% which renames the fields with the old name to the new name. Fields that
% are specified but not present will be silently ignored.
%
% See also COPYFIELDS, KEEPFIELDS, REMOVEFIELDS

% Copyright (C) 2014, Robert Oostenveld
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

if isempty(a)
  % this prevents problems if a is an empty double, i.e. []
  return
end

% these should be cell-arrays
if ischar(old)
  old = {old};
end
if ischar(new)
  new = {new};
end

if length(old)~=length(new)
  ft_error('the number of field names does not match between old and new');
end

% keep the fields that were not mentioned
b = keepfields(a, setdiff(fieldnames(a), old));
% copy the fields over with their new name
for i=1:length(old)
  if isfield(a, old{i})
    for j=1:numel(b)
      b(j).(new{i}) = a(j).(old{i});
    end
  end
end
