function SP=nifti2matlab(blob)
% NIFTI2MATLAB parses a NIFTI (currently NIFTI-1 only) header and generates a 
% corresponding MATLAB data structure. Use as
%
% NH = nifti2matlab(blob)
%
% where 'blob' needs to be of type uint8 or char. 
%
% Copyright (C) 2010, Stefan Klanke

warning('Trying to compile MEX file')
oldDir = pwd;
[srcDir, srcName] = fileparts(mfilename('fullpath'));
try
  cd(srcDir);			% we're now in private
  cd('../siemens');
  mex -outdir ../private nifti2matlab.c -I.
  cd(oldDir);
  SP = nifti2matlab(blob);
catch
  cd(oldDir);
  rethrow(lasterror)
end

