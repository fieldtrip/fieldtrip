function ft_test(varargin)

% FT_TEST performs selected FieldTrip test scripts or reports on previous test
% results from the dashboard.
%
% Use as
%   ft_test run             ...
%   ft_test moxunit_run     ...
%   ft_test report          ...
%   ft_test compare         ...
%
% ========= Running tests =========
%
% To execute a test and submit the results to the database, you would do
%   ft_test run
% to run all, or for a specific one
%   ft_test run test_bug46
%
% Test functions should not require any input arguments. Any output arguments will
% not be considered.
%
% Additional optional arguments are specified as key-value pairs and can include
%   dependency       = string
%   maxmem           = number (in bytes) or string such as 10GB
%   maxwalltime      = number (in seconds) or string such as HH:MM:SS
%   upload           = string, can be 'yes' or 'no' (default = 'yes')
%
% ========= Running MOxUnit tests =========
%
% To execute tests using MOxUNit, you would do
%   ft_test moxunit_run
%
% This feature is currently experimental, but should support the same 
% options as ft_test run (see above), and in addition:
%   xmloutput         = string, filename for JUnit-like XML file with test
%                       results (used for shippable CI)
%
% ========= Reporting on tests =========
%
% To print a table with the results on screen, you would do
%   ft_test result
% to show all, or for a specific one
%   ft_test result test_bug46
%
% Additional query arguments are specified as key-value pairs and can include
%   matlabversion    = string, for example 2016b
%   fieldtripversion = string
%   branch           = string
%   arch             = string, can be glnxa64, maci64. win32 or win64 
%   hostname         = string
%   user             = string
%
% ========= Comparing tests =========
%
% To print a table comparing different test results, you would do
%   ft_test compare matlabversion     ea3c2b9  314d186
%   ft_test compare fieldtripversion  2015b    2016b
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
% See also FT_VERSION

% Copyright (C) 2016-2017, Robert oostenveld
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

switch (varargin{1})
  case 'run'
    ft_test_run(varargin{:});
  case 'moxunit_run'
    ft_test_moxunit_run(varargin{:});
  case 'report'
    ft_test_report(varargin{:});
  case 'compare'
    ft_test_compare(varargin{:});
  otherwise
    error('unsupported command "%s"', varargin{1})
end % switch


