function y = struct(x, varargin)

% STRUCT Convert a config object into a structure object.

% Copyright (C) 2012-2015, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

if nargin==1
  % convert the config object into a regular MATLAB structure
  for i=1:numel(x)
    y(i) = struct(x(i).value);
  end
  y = reshape(y, size(x));
  % recurse into the structure to convert sub-configs into sub-structures
  key = fieldnames(y);
  for j=1:numel(y)
    for i=1:length(key)
      val = y(j).(key{i});
      if isa(val, 'config')
        y(j) = setfield(y(j), key{i}, struct(val));
      end
    end
  end
else
  % mimic the behaviour of the builtin MATLAB struct function
  if mod(nargin,2)
    error('Incorrect number of input arguments (should be key-value pairs)')
  end
  varargin = {x varargin{:}};
  key = varargin(1:2:end);
  val = varargin(2:2:end);

  y = struct();
  for i=1:length(key)
    y = setfield(y, key{i}, val{i});;
  end
end
