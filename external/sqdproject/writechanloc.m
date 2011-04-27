function writechanloc(sqdfile,locfile,format,channum)
% WRITECHANLOC(SQD_FILE, LOCATION_FILE, FORMAT, CHANNEL_LIST); Writes the location of 
% the sensors of the channels from the CHANNEL_LIST in the SQD_FILE to LOC_FILE 
% following the FORMAT structure
% where FORMAT can be:
%     'xyz' : (Default) Matlab/EEGLAB Cartesian coordinates (NOT EGI Cartesian).
%             x is toward nose; y is toward left ear; z is toward vertex
% Default values:
%   LOCATION_FILE = './testsqd.xyz'
%   FORMAT = 'xyz'
%   CHANNEL_LIST = 0:156
%  
% This is useful for analysing the data in EEGLAB toolbox
% Please see  http://sccn.ucsd.edu/eeglab for more details of the toolbox.
%
% Examples:
% If the path for the sqd-toolbox has been set properly and all of the files present 
% (see README), there should be a sample raw sqd-file "testsqd.sqd" for testing. The
% sqd-file "testsqd.sqd" has data for 192 channels, 10 seconds of data and sampled
% at 1000Hz. 
% 1. To read in channel location from "testsqd.sqd":
%   locationfile = 'testlocation.xyz';  % Filename for the location file to be created
%   sqdfilename  = 'testsqd.sqd';       % sqd-file to read the locations from
%   writechanloc(sqdfilename,locationfile); % Write location file
%   type(locationfile);                 % To verify
% 2. To read in channel location from "testsqd.sqd" for channels N1,N2,N3:
%   locationfile = 'testlocation.xyz';  % Filename for the location file to be created
%   sqdfilename  = 'testsqd.sqd';       % sqd-file to read the locations from
%   N1 = 1;N2 = 100;N3 = 101;           % Init channels to read
%   writechanloc(sqdfilename,locationfile,'xyz',[N1,N2,N3]); % Write location file
%   type(locationfile);                 % To verify
%
% See also,
% @chanhandle/writechanloc

% Note:
% This function is just a wrapper around @chanhandle/writechanloc
% Error correction is taken care by respective functions

% Init input arguments
switch nargin
case 0
    error('Atleast one input - SQDFILENAME required');
case 1
    locfile = 'testsqd.xyz';
    format = 'xyz';
    channum = 0:156;
case 2
    format = 'xyz';
    channum = 0:156;
case 3
    channum = 0:156;
end;

chant = chanhandle(sqdfile,channum);    % get chanhandle
writechanloc(chant,locfile,format);     % Write file
