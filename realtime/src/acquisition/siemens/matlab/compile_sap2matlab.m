% SAP2MATLAB parses Siemens ASCII protocol data and generates a 
% corresponding MATLAB data structure.
%
% This function is used for de-serialising the header information
% from a FieldTrip buffer containing fMRI data from a Siemens scanner.

% Copyright (C) 2010, Stefan Klanke,
% 	Modified by Tim van Mourik, 2014
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

warning('Trying to compile MEX file')
oldDir = pwd;
try
    %This is where the header files are located
    options = '-I../include';
    %These two C-files are copiled 
    fileNames = '../src/sap2matlab.c ../src/siemensap.c ';
    %the compilation does not depend on external libraries
    libraries = [];
    eval(['mex ' fileNames, libraries, options]);
catch
    rethrow(lasterror)
end

%% Test
%  This file contains an example string that will be parsed by sap2matlab
load('mrprotString.mat', 'apstr');
S = sap2matlab(apstr);


%%
file = '../Projects/TestReadMrprot/mrprot_triotim.txt';
f = fopen(file, 'r');
textData = fread(f);
fclose(f);

S = sap2matlab(sprintf('%s', textData));


