function status = webwrite(varargin)

% WEBWRITE is a drop-in replacement for the function with the same
% name that was introduced in MATLAB 2014b. This function is only
% partially compatible with the original.
%
% This requires that curl is available on the command-line.

% Copyright (C) 2017, Robert Oostenveld
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

url = varargin{1};

if nargin>2
  data = varargin{2};
else
  data = '';
end

if nargin>2
  options = varargin{3};
else
  options = {};
end

if isstruct(data)
  % convert it to JSON
  data = struct2json(data);
  data(data==10) = []; % remove newlines
  data(data==32) = []; % remove spaces
  cmd = sprintf('curl -H "Content-Type: application/json; charset=UTF-8" -X POST -d ''%s'' %s', data, url);
elseif ischar(data)
  % send it as it is
  cmd = sprintf('curl -H "Content-Type: application/text; charset=UTF-8" -X POST -d ''%s'' %s', data, url);
else
  error('unsupported data');
end

[status, str] = system(cmd);
