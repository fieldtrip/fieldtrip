function [val, remaining] = keyval(key, varargin)

% KEYVAL returns the value that corresponds to the requested key in a
% key-value pair list of variable input arguments
%
% Use as
%   [val] = keyval(key, varargin)
%
% See also VARARGIN

% Undocumented option
%   [val] = keyval(key, varargin, default)

% Copyright (C) 2005-2007, Robert Oostenveld
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

% what to return if the key is not found
emptyval = [];

if nargin==3 && iscell(varargin{1})
  emptyval = varargin{2};
  varargin = varargin{1};
end

if nargin==2 && iscell(varargin{1})
  varargin = varargin{1};
end

if mod(length(varargin),2)
  ft_error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
end

% the 1st, 3rd, etc. contain the keys, the 2nd, 4th, etc. contain the values
keys = varargin(1:2:end);
vals = varargin(2:2:end);

% the following may be faster than cellfun(@ischar, keys)
valid = false(size(keys));
for i=1:numel(keys)
  valid(i) = ischar(keys{i});
end

if ~all(valid)
  ft_error('optional input arguments should come in key-value pairs, the optional input argument %d is invalid (should be a string)', i);
end

hit = find(strcmpi(key, keys));
if isempty(hit)
  % the requested key was not found
  val = emptyval;
elseif length(hit)==1  
  % the requested key was found
  val = vals{hit};
else
  ft_error('multiple input arguments with the same name');
end

if nargout>1
  % return the remaining input arguments with the key-value pair removed
  keys(hit) = [];
  vals(hit) = [];
  remaining = cat(1, keys(:)', vals(:)');
  remaining = remaining(:)';
end
