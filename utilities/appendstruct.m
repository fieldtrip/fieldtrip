function s = appendstruct(varargin)

% APPENDSTRUCT appends a structure or a struct-array to another structure or
% struct-array. It also works if the initial structure is an empty structure or an
% empty double array. It also works if the input structures have different fields.
%
% Use as
%   ab = appendstruct(a, b)
%
% See also PRINTSTRUCT, MERGESTRUCT, COPYFIELDS, KEEPFIELDS, REMOVEFIELDS

% Copyright (C) 2015-2022, Robert Oostenveld
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

narginchk(2,inf);

for i=1:nargin
  assert(isstruct(varargin{i}) || isempty(varargin{i}), 'input argument %d should be empty or a structure', i);
end

if nargin>2
  % use recursion to append multiple structures
  s = varargin{1};
  for i=2:nargin
    s = appendstruct(s, varargin{i});
  end
  return
else
  % continue with the code below to append two structures
  s1 = varargin{1};
  s2 = varargin{2};
end

if isempty(s1) && isempty(s2)
  % this results in a 0x0 empty struct array with no fields
  s = struct([]);
elseif isempty(s1) && ~isempty(s2)
  % return only the second one
  s = s2;
elseif isempty(s2) && ~isempty(s1)
  % return only the first one
  s = s1;
else
  % concatenate the second structure to the first
  fn1 = fieldnames(s1);
  fn2 = fieldnames(s2);
  % find the fields that are missing in either one
  missing1 = setdiff(union(fn1, fn2), fn1);
  missing2 = setdiff(union(fn1, fn2), fn2);
  % add the missing fields
  for i=1:numel(missing1)
    s1(1).(missing1{i}) = [];
  end
  for i=1:numel(missing2)
    s2(1).(missing2{i}) = [];
  end
  s = cat(1, s1(:), s2(:));
end
