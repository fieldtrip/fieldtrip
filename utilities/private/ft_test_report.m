function result = ft_test_report(varargin)

% FT_TEST_REPORT

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

narginchk(1, inf);
command = varargin{1};
assert(isequal(command, 'report'));
varargin = varargin(2:end);

optbeg = find(ismember(varargin, {'matlabversion', 'fieldtripversion', 'user', 'hostname', 'branch', 'arch', 'showdate', 'showid', 'timeout'}));
if ~isempty(optbeg)
  optarg   = varargin(optbeg:end);
  varargin = varargin(1:optbeg-1);
else
  optarg = {};
end

% varargin contains the file (or files) to test
% optarg contains the command-specific options
showdate = ft_getopt(optarg, 'showdate', false);  if ischar(showdate), showdate = istrue(showdate);    end
showid   = ft_getopt(optarg, 'showid', false);    if ischar(showid),   showid   = istrue(showid);      end
timeout  = ft_getopt(optarg, 'timeout', 30);      if ischar(timeout),  timeout  = str2double(timeout); end

% construct the query string that will be passed in the URL
query = '?';
queryparam = {'matlabversion', 'fieldtripversion', 'hostname', 'user', 'branch', 'arch'};
for i=1:numel(queryparam)
  val = ft_getopt(optarg, queryparam{i});
  if ~isempty(val)
    query = [query sprintf('%s=%s&', queryparam{i}, val)];
  end
end

options = weboptions('ContentType', 'json', 'Timeout', timeout); % this returns the result as MATLAB structure
url = 'http://dashboard.fieldtriptoolbox.org/api/';

if isempty(varargin)
  result = webread([url query], options);
  assert(~isempty(result), 'no results were returned');
  result = mergecellstruct(result);
else
  results = cell(size(varargin));
  for i=1:numel(varargin)
    result = webread([url query sprintf('&functionname=%s', varargin{i})], options);
    assert(~isempty(result), 'no results were returned for %s', varargin{i});
    results{i} = mergecellstruct(result);
  end
  % merge all results
  result = mergecellstruct(results);
end

% rename the automatically added fields
result = renamefields(result, 'x_id', 'id');
result = renamefields(result, 'createDate', 'date');

% remove some of the fields
if ~showid
  result = removefields(result, 'id');
end
if ~showdate
  result = removefields(result, 'date');
end

% convert the struct-array to a table
table = struct2table(result);
fprintf('%s\n', table{:});

