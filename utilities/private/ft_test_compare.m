function [summary] = ft_test_compare(varargin)

% FT_TEST_COMPARE

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

narginchk(2, inf);
command = varargin{1};
feature = varargin{2};
assert(isequal(command, 'compare'));
varargin = varargin(3:end);

optbeg = find(ismember(varargin, {'matlabversion', 'fieldtripversion', 'user', 'hostname', 'branch', 'arch'}));
if ~isempty(optbeg)
  optarg   = varargin(optbeg:end);
  varargin = varargin(1:optbeg-1);
else
  optarg = {};
end

% varargin contains the file (or files) to test
% optarg contains the command-specific options

% construct the query string that will be passed in the URL
query = '?';
queryparam = {'matlabversion', 'fieldtripversion', 'hostname', 'user', 'branch', 'arch'};
for i=1:numel(queryparam)
  val = ft_getopt(optarg, queryparam{i});
  if ~isempty(val)
    query = [query sprintf('%s=%s&', queryparam{i}, val)];
  end
end

options = weboptions('ContentType','json'); % this returns the results as MATLAB structure
url = 'http://dashboard.fieldtriptoolbox.org/api/';

results = cell(size(varargin));
functionname = {};
for i=1:numel(varargin)
  result = webread([url query sprintf('&%s=%s', feature, varargin{i})], options);
  
  % the documents in the mongoDB database might not fully consistent, in which case they are returned as cell-array containing different structures
  % merge all stuctures into a single struct-array
  result = mergecellstruct(result);
  
  assert(~isempty(result), 'no results were returned for %s %s', feature, varargin{i});
  functionname = cat(1, functionname(:), {result.functionname}');
  results{i} = result;
end

% find the joint set of all functions
functionname = unique(functionname);

% represent a summary of all results in a struct-array
summary = struct();
for i=1:numel(functionname)
  summary(i).function = functionname{i};
  for j=1:numel(varargin)
    sel = find(strcmp({results{j}.functionname}, functionname{i}));
    fn = fixname(varargin{j}, 'X_base64encode_X');
    summary(i).(fn) = haspassed(results{j}, sel);
  end % for each functionname
end % for each of the features

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to convert the boolean result into a string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = haspassed(result, index)
if isempty(index)
  str = [];
elseif isfield(result, 'passed') && all(istrue([result(index).passed]))
  str = 'passed';
elseif isfield(result, 'passed') && all(~istrue([result(index).passed]))
  str = 'failed';
elseif isfield(result, 'result') && all(istrue([result(index).result]))
  str = 'passed';
elseif isfield(result, 'result') && all(~istrue([result(index).result]))
  str = 'failed';
else
  % there appear to be multiple representations of the same test, but they are not consistent
  str = 'ambiguous';
end

