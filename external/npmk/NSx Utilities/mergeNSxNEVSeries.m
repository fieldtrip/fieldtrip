%
%   mergeNSxNEV()
%
%   This function loads two NSx and NEV files and it will combine them
%   together into one file. The resulting file will be saved as new NSx
%   and NEV files onto the disk.  To combine two NSx and NEV files into 
%   indivual NSx and NEV variables in MATLAB see combineNSxNEV. The time
%   difference between the two sets of recordings is removed. To determine
%   the time differnce between the two data files, use
%   NSx.MetaTags.DateTimeRaw or NEV.MetaTags.DateTimeRaw variables.
%
%
%   filename1:  The name of the first NSx file. This input is optional. In
%               its absense, a dialog will open and will prompt the user to
%               select an NSx file.
%               (OPTIONAL)
%
%   filename2:  The name of the second NSx file. This input is also
%               optional. In its absense, a dialog will open and will
%               prompt the user to select an NSx file.
%               (OPTIONAL)
%
%   Example: 
%   
%   mergeNSxNEV('c:\data\saveddata1.ns5', 'c:\data\saveddata2.ns5');
%
%   The above example reads the two files (full path needed)
%   c:\data\saveddata1.ns5 and c:\data\saveddata2.ns5 and their corresponding
%   NEV files (saveddata1.nev and saveddata2.nev) in the same folder and
%   merges them into a single file called firstrecording001-combined.ns2
%   and firstrecording001-combined.nev.
%
%   Kian Torab
%   kian@blackrockmicro.com
%   Blackrock Microsystems
%   Version 1.0.0.0
%   March 31, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mergeNSxNEV(varargin)

if nargin >= 1
    disp('Too many input arguments.');
    return;
end

%% Opening the file
% Getting the name of the first file from the user
[gfFileName gfPathName] = getFile;
[filePath, fileName, fileExt] = fileparts(gfFileName);


