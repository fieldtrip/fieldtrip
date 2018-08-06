function val = ft_getopt(opt, key, default, emptymeaningful)

% FT_GETOPT gets the value of a specified option from a configuration structure
% or from a cell-array with key-value pairs.
%
% Use as
%   val = ft_getopt(s, key, default, emptymeaningful)
% where the input values are
%   s               = structure or cell-array
%   key             = string
%   default         = any valid MATLAB data type (optional, default = [])
%   emptymeaningful = boolean value (optional, default = 0)
%
% If the key is present as field in the structure, or as key-value pair in the
% cell-array, the corresponding value will be returned.
%
% If the key is not present, ft_getopt will return the default, or an empty array
% when no default was specified.
%
% If the key is present but has an empty value, then the emptymeaningful flag
% specifies whether the empty value or the default value should be returned.
% If emptymeaningful==true, then the empty array will be returned.
% If emptymeaningful==false, then the specified default will be returned.
%
% See also FT_SETOPT, FT_CHECKOPT

% Copyright (C) 2011-2012, Robert Oostenveld
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

if nargin<3
  default = [];
end

if nargin < 4
  emptymeaningful = 0;
end

if isa(opt, 'struct') || isa(opt, 'config')
  % get the key-value from the structure
  fn = fieldnames(opt);
  if ~any(strcmp(key, fn))
    val = default;
  else
    val = opt.(key);
  end
  
elseif isa(opt, 'cell')
  % get the key-value from the cell-array
  if mod(length(opt),2)
    error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
  end
  
  % the 1st, 3rd, etc. contain the keys, the 2nd, 4th, etc. contain the values
  keys = opt(1:2:end);
  vals = opt(2:2:end);
  
  % the following may be faster than cellfun(@ischar, keys)
  valid = false(size(keys));
  for i=1:numel(keys)
    valid(i) = ischar(keys{i});
  end
  
  if ~all(valid)
    error('optional input arguments should come in key-value pairs, the optional input argument %d is invalid (should be a string)', i);
  end
  
  hit = find(strcmpi(key, keys));
  if isempty(hit)
    % the requested key was not found
    val = default;
  elseif length(hit)==1
    % the requested key was found
    val = vals{hit};
  else
    error('multiple input arguments with the same name');
  end
  
elseif isempty(opt)
  % no options are specified, return default
  val = default;
end % isstruct or iscell or isempty

if isempty(val) && ~isempty(default) && ~emptymeaningful
  % use the default value instead of the empty input that was specified:
  % this applies for example if you do functionname('key', []), where
  % the empty is meant to indicate that the user does not know or care
  % what the value is
  val = default;
end
