function results = ft_test_result(varargin)

% FT_TEST_RESULT checks the status of selected test scripts on the FieldTrip dashboard
%
% Use as
%   ft_test_result
%
% Query arguments are specified as key-value pairs and can include
%   matlabversion    = string
%   fieldtripversion = string
%   hostname         = string
%   user             = string
%
% To get a list of all distinct values of a certain parameter, you specify
%   distinct         = string
%
% See also FT_TEST_RUN, FT_VERSION

% Copyright (C) 2016, Robert oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/donders/fieldtrip
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

% set the default
command = 'query';

if nargin>0 && isequal(varargin{1}, 'compare')
  % compare rev1 rev2
  command   = varargin{1};
  revision1 = varargin{2};
  revision2 = varargin{3};
  varargin = varargin(4:end);
end

% construct the query string
query = '?';

queryparam = {'matlabversion', 'fieldtripversion', 'hostname', 'user', 'functionname', 'distinct'};
for i=1:numel(queryparam)
  val = ft_getopt(varargin, queryparam{i});
  if ~isempty(val)
    query = [query sprintf('%s=%s&', queryparam{i}, val)];
  end
end

options = weboptions('ContentType','json'); % this returns the results as MATLAB structure

switch command
  case 'query'
    results = webread(['http://dashboard.fieldtriptoolbox.org/api/' query], options);

  case 'compare'
    dashboard1 = webread(['http://dashboard.fieldtriptoolbox.org/api/' query sprintf('&fieldtripversion=%s', revision1)], options);
    dashboard2 = webread(['http://dashboard.fieldtriptoolbox.org/api/' query sprintf('&fieldtripversion=%s', revision2)], options);
    assert(~isempty(dashboard1), 'no tests were returned for the first revision');
    assert(~isempty(dashboard2), 'no tests were returned for the second revision');
    functionname1 = {dashboard1.functionname};
    functionname2 = {dashboard2.functionname};
    functionname = unique(cat(2, functionname1, functionname2));
    n = max(cellfun(@length, functionname));
    for i=1:numel(functionname)
      sel1 = find(strcmp(functionname1, functionname{i}));
      sel2 = find(strcmp(functionname2, functionname{i}));
      res1 = getresult(dashboard1, sel1);
      res2 = getresult(dashboard2, sel2);
      fprintf('%s : %s in %s, %s in %s\n', padto(functionname{i}, n), res1, revision1, res2, revision2);
    end % for
    
  otherwise
    error('unsupported command "%s"', command);
end % switch command

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = padto(str, n)
if n>length(str)
  str = [str repmat(' ', [1 n-length(str)])];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = getresult(dashboard, sel)
if isempty(sel)
  str = 'missing';
elseif all(istrue([dashboard(sel).result]))
  str = 'passed';
elseif all(~istrue([dashboard(sel).result]))
  str = 'failed';
else
  str = 'ambiguous';
end