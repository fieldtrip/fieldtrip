function opt = ft_deleteopt(opt, keys)

% FT_DELETEOPT removes one or more keys and their corresponding values from a
% configuration structure or from a cell-array with key-value pairs. It will 
% ignore keys that are not in the structure/cell-array.
%
% Use as
%   s = ft_deleteopt(s, keys)
% where s is a structure or a cell-array, and keys is a cell-array or char.
%
% See also FT_GETOPT, FT_CHECKOPT, FT_SETOPT

% Copyright (C) 2025, JAN MATHIJS SCHOFFELEN
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

if isa(keys, 'char')
  keys = {keys};
end

if isa(opt, 'struct')
  % just remove the fields from the struct
  opt = removefields(opt, keys);
elseif isa(opt, 'cell')
  sel = true(1, numel(opt));
  for i=1:numel(keys)
    thiskey = find(strcmp(opt(1:2:end), keys{i}));
    if ~isempty(thiskey)
      sel(2*(thiskey-1) + 1)   = false; % deselect the key
      sel(2*(thiskey-1) + 2) = false; % deselect the value
    end
  end
  opt = opt(sel);
end % isstruct or iscell
