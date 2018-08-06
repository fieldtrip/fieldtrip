function [s] = removefields(s, fields, varargin)

% REMOVEFIELDS makes a selection of the fields in a structure
%
% Use as
%   s = removefields(s, fields);
%
% See also KEEPFIELDS, COPYFIELDS

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

if isempty(s)
  % this prevents problems if s is an empty double, i.e. []
  return
end

% get the optional arguments
recursive = ft_getopt(varargin, 'recursive', false);

if ischar(fields)
  fields = {fields};
elseif ~iscell(fields)
  ft_error('fields input argument must be a cell array of strings or a single string');
end

remove = intersect(fieldnames(s), fields);
for i=1:numel(remove)
  s = rmfield(s, remove{i});
end

if recursive
  fn = fieldnames(s);
  for i=1:numel(fn)
    if isstruct(s.(fn{i}))
      s.(fn{i}) = removefields(s.(fn{i}), fields, varargin{:});
    end
  end
end
