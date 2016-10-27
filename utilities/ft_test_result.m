function results = ft_test_result(varargin)

% FT_TEST_RESULT checks the status of selected test scripts on the FieldTrip dashboard
%
% To get all dashboard results as a structure array, you would do
%   result = ft_test_result
%
% To print a table with the results on screen, you would do
%   ft_test_result comparerevision ea3c2b9 314d186
%   ft_test_result comparematlab   2015b   2016b
%   ft_test_result comparemexext   mexw32  mexw64
%   ft_test_result compareos       osx     windows
%
% Additional query arguments are specified as key-value pairs and can include
%   matlabversion    = string
%   fieldtripversion = string
%   hostname         = string
%   user             = string
%   branch           = string
%   result           = string
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

if nargin>0
  if isequal(varargin{1}, 'compare') || isequal(varargin{1}, 'comparerevision')
    % comparerevision rev1 rev2
    command  = 'comparerevision';
    arg1     = varargin{2};
    arg2     = varargin{3};
    varargin = varargin(4:end);
  elseif isequal(varargin{1}, 'matlab') || isequal(varargin{1}, 'comparematlab')
    % comparematlab ver1 ver2
    command  = 'comparematlab';
    arg1     = varargin{2};
    arg2     = varargin{3};
    varargin = varargin(4:end);
  end
end

% construct the query string
query = '?';

% the 'distict' parameter is mutually exclusive with all others
queryparam = {'matlabversion', 'fieldtripversion', 'hostname', 'user', 'functionname', 'result', 'branch', 'distinct'};
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
    
  case 'comparerevision'
    dashboard1 = webread(['http://dashboard.fieldtriptoolbox.org/api/' query sprintf('&fieldtripversion=%s', arg1)], options);
    dashboard2 = webread(['http://dashboard.fieldtriptoolbox.org/api/' query sprintf('&fieldtripversion=%s', arg2)], options);
    assert(~isempty(dashboard1), 'no tests were returned for the first revision');
    assert(~isempty(dashboard2), 'no tests were returned for the second revision');
    functionname1 = {dashboard1.functionname};
    functionname2 = {dashboard2.functionname};
    functionname = unique(cat(2, functionname1, functionname2));
    n = max(cellfun(@length, functionname));
    line = cell(size(functionname));
    order = nan(size(functionname));
    for i=1:numel(functionname)
      sel1 = find(strcmp(functionname1, functionname{i}));
      sel2 = find(strcmp(functionname2, functionname{i}));
      res1 = getresult(dashboard1, sel1);
      res2 = getresult(dashboard2, sel2);
      line{i} = sprintf('%s : %s in %s, %s in %s\n', padto(functionname{i}, n), res1, arg1, res2, arg2);
      
      % determine the order
      if     strcmp(res1, 'passed') && strcmp(res2, 'passed')
        order(i) = 1;
      elseif strcmp(res1, 'failed') && strcmp(res2, 'failed')
        order(i) = 2;
      elseif strcmp(res1, 'missing') && strcmp(res2, 'passed')
        order(i) = 3;
      elseif strcmp(res1, 'passed') && strcmp(res2, 'missing')
        order(i) = 4;
      elseif strcmp(res1, 'missing') && strcmp(res2, 'failed')
        order(i) = 5;
      elseif strcmp(res1, 'failed') && strcmp(res2, 'missing')
        order(i) = 6;
      elseif strcmp(res1, 'failed') && strcmp(res2, 'passed')
        order(i) = 7;
      elseif strcmp(res1, 'passed') && strcmp(res2, 'failed')
        order(i) = 8;
      else
        order(i) = 9;
      end
    end % for each functionname
    
    [order, index] = sort(order);
    line = line(index);
    for i=1:length(line)
      fprintf(line{i});
    end
    
  case 'comparematlab'
    dashboard1 = webread(['http://dashboard.fieldtriptoolbox.org/api/' query sprintf('&matlabversion=%s', arg1)], options);
    dashboard2 = webread(['http://dashboard.fieldtriptoolbox.org/api/' query sprintf('&matlabversion=%s', arg2)], options);
    assert(~isempty(dashboard1), 'no tests were returned for the first matab version');
    assert(~isempty(dashboard2), 'no tests were returned for the second matab version');
    functionname1 = {dashboard1.functionname};
    functionname2 = {dashboard2.functionname};
    functionname = unique(cat(2, functionname1, functionname2));
    n = max(cellfun(@length, functionname));
    line = cell(size(functionname));
    order = nan(size(functionname));
    for i=1:numel(functionname)
      sel1 = find(strcmp(functionname1, functionname{i}));
      sel2 = find(strcmp(functionname2, functionname{i}));
      res1 = getresult(dashboard1, sel1);
      res2 = getresult(dashboard2, sel2);
      line{i} = sprintf('%s : %s in %s, %s in %s\n', padto(functionname{i}, n), res1, arg1, res2, arg2);
      
      % determine the order
      if     strcmp(res1, 'passed') && strcmp(res2, 'passed')
        order(i) = 1;
      elseif strcmp(res1, 'failed') && strcmp(res2, 'failed')
        order(i) = 2;
      elseif strcmp(res1, 'missing') && strcmp(res2, 'passed')
        order(i) = 3;
      elseif strcmp(res1, 'passed') && strcmp(res2, 'missing')
        order(i) = 4;
      elseif strcmp(res1, 'missing') && strcmp(res2, 'failed')
        order(i) = 5;
      elseif strcmp(res1, 'failed') && strcmp(res2, 'missing')
        order(i) = 6;
      elseif strcmp(res1, 'failed') && strcmp(res2, 'passed')
        order(i) = 7;
      elseif strcmp(res1, 'passed') && strcmp(res2, 'failed')
        order(i) = 8;
      else
        order(i) = 9;
      end
    end % for each functionname
    
    [order, index] = sort(order);
    line = line(index);
    for i=1:length(line)
      fprintf(line{i});
    end
    
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
