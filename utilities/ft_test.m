function status = ft_test(varargin)

% FT_TEST performs selected FieldTrip test scripts or reports on previous test
% results from the dashboard.
%
% Use as
%   ft_test run     ...
%   ft_test report  ...
%   ft_test compare ...
%
% === Running tests ===
%
% To execute a test and submit the results to the database, you would do
%   ft_test run test_bug46
%
% Test functions should not require any input arguments. Any output arguments will
% not be considered.
%
% Additional optional arguments are specified as key-value pairs and can include
%   dependency       = string
%   maxmem           = number (in bytes) or string
%   maxwalltime      = number (in seconds) string
%
% === Reporting on tests ===
%
% To print a table with the results on screen, you would do
%   ft_test result
%
% Additional query arguments are specified as key-value pairs and can include
%   matlabversion    = string
%   fieldtripversion = string
%   hostname         = string
%   user             = string
%   branch           = string
%
% === Comparing tests ===
%
% To print a table with the results on screen, you would do
%   ft_test comparerevision ea3c2b9 314d186
%   ft_test comparematlab   2015b   2016b
%   ft_test comparemexext   mexw32  mexw64
%   ft_test compareos       osx     windows
%
% Additional query arguments are specified as key-value pairs and can include
%   matlabversion    = string
%   fieldtripversion = string
%   hostname         = string
%   user             = string
%   branch           = string
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
    ft_test_run(varargin{1:end});
  case 'report'
    ft_test_report(varargin{1:end});
  case {'comparerevision', 'comparematlab'}
    ft_test_compare(varargin{1:end});
  otherwise
    error('unsupported command "%s"', varargin{1})
end % switch


