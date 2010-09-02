function SP=sap2matlab(blob)
% SAP2MATLAB parses Siemens ASCII protocol data and generates a 
% corresponding MATLAB data structure. Use as
%
% SP = sap2matlab(blob)
%
% where 'blob' needs to be of type uint8 or char. This function
% is currently used for de-serialising the header information
% from a FieldTrip buffer containing fMRI data from a Siemens
% scanner.
%
% The sources for the MEX-file are in the realtime/siemens directory,
% but calling the .m-file will try to compile automatically.
%
% Copyright (C) 2010, Stefan Klanke

warning('Trying to compile MEX file')
oldDir = pwd;
[srcDir, srcName] = fileparts(mfilename('fullpath'));
try
  cd(srcDir);			% we're now in private
  cd('../../realtime/datasource/siemens');
  mex -outdir ../../fileio/private sap2matlab.c siemensap.c -I.
  cd(oldDir);
  SP = sap2matlab(blob);
catch
  cd(oldDir);
  rethrow(lasterror)
end

