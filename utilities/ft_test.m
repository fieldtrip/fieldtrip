function [result] = ft_test(varargin)

% FT_TEST performs selected FieldTrip test scripts, finds and updates the dependencies of test scripts, finds which high-level FieldTrip functions are not tested, or reports on previous test
% results from the dashboard database.
%
% Use as
%   ft_test inventorize        ...
%   ft_test run                ...
%   ft_test find_dependency    ...
%   ft_test update_dependency  ...
%   ft_test untested_functions ...
%   ft_test moxunit_run        ... % this is obsolete
%   ft_test report             ... % this is obsolete
%   ft_test compare            ... % this is obsolete
% 
% ========= INVENTORIZE =========
% 
% To list the tests based on their dependencies, you would do 
%   ft_test inventorize
% to list all test functions, or
%   ft_test inventorize data no
% to list test functions that don't need any external data to run.
%  
% Additional optional arguments are specified as key-value pairs and can include
%   dependency       = string or cell-array of strings
%   upload           = string, upload test results to the dashboard, can be 'yes' or 'no' (default = 'yes')
%   dccnpath         = string, allow files to be read from the DCCN path, can be 'yes' or 'no' (default is automatic)
%   maxmem           = number (in bytes) or string such as 10GB
%   maxwalltime      = number (in seconds) or string such as HH:MM:SS
%   data             = string or cell-array of strings with 'no', 'public' and/or 'private'
%   sort             = string, can be 'alphabetical', 'walltime', 'mem' or 'random' (default = 'alphabetical')
%   returnerror      = string, whether give an error upon detecting a failed script, can be 'immediate', 'final', 'no' (default = 'no')
% 
% ========= RUN =========
%
% To execute a test and submit the results to the dashboard database, you would do
%   ft_test run
% to run all test functions, or
%   ft_test run test_bug46
% to run a selected test.
%
% Test functions should not require any input arguments. Any output arguments will
% not be considered.
%
% Additional optional arguments are specified as key-value pairs and can include
%   dependency       = string or cell-array of strings
%   upload           = string, upload test results to the dashboard, can be 'yes' or 'no' (default = 'yes')
%   dccnpath         = string, allow files to be read from the DCCN path, can be 'yes' or 'no' (default is automatic)
%   maxmem           = number (in bytes) or string such as 10GB
%   maxwalltime      = number (in seconds) or string such as HH:MM:SS
%   data             = string or cell-array of strings with 'no', 'public' and/or 'private'
%   sort             = string, can be 'alphabetical', 'walltime', 'mem' or 'random' (default = 'alphabetical')
%   returnerror      = string, whether give an error upon detecting a failed script, can be 'immediate', 'final', 'no' (default = 'no')
% 
% ========= FIND_DEPENDENCY =========
% 
% To find on what functions test scripts depend on, you would do
%   ft_test find_dependency test_bug46
% to find on what functions test_bug46 depends on.
% 
% It outputs:
%   inlist   = Nx1 cell-array, describes the rows and lists the test scripts
%   outlist  = 1xM cell-array, describes the columns and lists the dependencies
%   depmat   = NxM dependency matrix, see below
% 
% The dependency matrix contains the following values:
%  - 0 if there is no dependency
%  - 2 for a direct dependency
% 
% ========= UPDATE_DEPENDENCY =========
% 
% To update the DEPENDENCY header in a specific test script, you would do:
%   ft_test update_dependency test_bug46
% 
% ========= UNTESTED_FUNCTIONS =========
% 
% To find FieldTrip high-level functions not tested by any test scripts,
% you would do
%   ft_test untested_functions
%
% ========= MOXUNIT_RUN =========
%
% To execute tests using MOxUNit, you would do
%   ft_test moxunit_run
%
% This feature is still experimental, but should support the same
% options as ft_test run (see above), and in addition:
%   xmloutput         = string, filename for JUnit-like XML file with test
%                       results (used for shippable CI).
%   exclude_if_prefix_equals_failed = string, if set to false (or 'no')
%                       then tests are also run if their filename starts
%                       with 'failed'. If set to true (or 'yes'), which is
%                       the default, then filenames starting with 'failed'
%                       are skipped.
%
% ========= REPORT =========
%
% To print a table with the results on screen, you would do
%   ft_test report
% to show all, or for a specific one
%   ft_test report test_bug46
%
% Additional query arguments are specified as key-value pairs and can include
%   matlabversion    = string, for example 2016b
%   fieldtripversion = string
%   branch           = string
%   arch             = string, can be glnxa64, maci64. win32 or win64
%   hostname         = string
%   user             = string
%
% Optionally, you may capture the output to get the results as a Matlab table
% array, in which case they are not automatically displayed.
%   rslt = ft_test('report', 'fieldtripversion', 'cef3396');
%
% ========= COMPARE =========
%
% To print a table comparing different test results, you would do
%   ft_test compare matlabversion     2015b    2016b
%   ft_test compare fieldtripversion  ea3c2b9  314d186
%   ft_test compare arch              glnxa64  win32
%
% Additional query arguments are specified as key-value pairs and can include
%   matlabversion    = string, for example 2016b
%   fieldtripversion = string
%   branch           = string
%   arch             = string, can be glnxa64, maci64. win32 or win64
%   hostname         = string
%   user             = string
%
% See also DCCNPATH, FT_VERSION

% Copyright (C) 2016-2024, Robert Oostenveld
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

% This function is designed to be executed like this
%   ft_test run test_bug46
% but you can also execute it like this
%   ft_test('run', 'test_bug46')
% which is required if you want to get the result as output argument.

if nargin==0
  help(mfilename)
  return
end

% ensure that all input arguments are strings, required for maxwalltime etc.
for i=1:numel(varargin)
  if isnumeric(varargin{i})
    varargin{i} = num2str(varargin{i});
  end
end

switch (varargin{1})
  case 'run'
    result = ft_test_run(varargin{:});
  case 'inventorize'
    result = ft_test_run(varargin{:}); % this uses the same code as 'run'
  case 'find_dependency'
    [outlist, depmat] = ft_test_find_dependency(varargin{:});
  case 'update_dependency'
    [outlist, depmat] = ft_test_find_dependency(varargin{:});
    ft_test_update_dependency(depmat, varargin(2:end), outlist);
  case 'untested_functions'
    ft_test_untested_functions(varargin{:});
  case 'report'
    result = ft_test_report(varargin{:});
  case 'compare'
    result = ft_test_compare(varargin{:});
  case 'moxunit_run'
    result = ft_test_moxunit_run(varargin{:});
  otherwise
    ft_error('unsupported command "%s"', varargin{1})
end % switch

switch (varargin{1})
  case 'find_dependency'
    fprintf('\n');
    fprintf('The test scripts: ');
    fprintf('<strong>%s</strong> ', varargin{2:end});
    fprintf('depend on: ');
    fprintf('<strong>%s</strong> ', outlist{:});
    fprintf('\n\n');

  case 'update_dependency'
    return
  
  case 'untested_functions'
    return

otherwise
    if ~nargout
      % show it on screen, do not return 'ans'
      if isempty(result)
        fprintf('Results are empty\n');
      else                  
         printstruct_as_table(result);                
      end
      clear result
    else
      % convert it to a proper MATLAB table
      result = struct2table(result);
    end
end
