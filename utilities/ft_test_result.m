function results = ft_test_result(varargin)

% FT_TEST_RESULT checks the status of selected test scripts on the FieldTrip dashboard
%
% Use as
%   ft_test_result
%
% Additional arguments are specified as key-value pairs and can include
%   matlabversion    = string
%   fieldtripversion = string
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

optbeg = find(ismember(varargin, {'matlabversion', 'fieldtripversion'}));
if ~isempty(optbeg)
  optarg = varargin(optbeg:end);
  varargin = varargin(1:optbeg-1);
else
  optarg = {};
end

% get the optional input arguments
matlabversion    = ft_getopt(optarg, 'matlabversion', {});
fieldtripversion = ft_getopt(optarg, 'fieldtripversion', inf);

if ischar(matlabversion)
  % this should be a cell-array
  matlabversion = {matlabversion};
end

if ischar(fieldtripversion)
  % this should be a cell-array
  fieldtripversion = {fieldtripversion};
end

options = weboptions('ContentType','json'); % this returns the results as MATLAB structure
results = webread('http://dashboard.fieldtriptoolbox.org/test', options);

