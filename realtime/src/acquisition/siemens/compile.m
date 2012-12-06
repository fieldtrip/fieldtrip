function SP=sap2matlab(blob)

% SAP2MATLAB parses Siemens ASCII protocol data and generates a 
% corresponding MATLAB data structure.
%
% Use as
%   SP = sap2matlab(blob)
%
% where 'blob' needs to be of type uint8 or char. This function
% is currently used for de-serialising the header information
% from a FieldTrip buffer containing fMRI data from a Siemens
% scanner.

% Copyright (C) 2010, Stefan Klanke
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

warning('Trying to compile MEX file')
oldDir = pwd;
[srcDir, srcName] = fileparts(mfilename('fullpath'));
try
  cd(srcDir);			% we're now in private
  cd('../siemens');
  mex sap2matlab.c siemensap.c -I.
  cd(oldDir);
  SP = sap2matlab(blob);
catch
  cd(oldDir);
  rethrow(lasterror)
end

