function separatePausedNSx(varargin)

% separatePausedNsx    Saves paused data segments from a single NSx file 
%                      as individual NSx files
%
%    separatePausedNSx(FILENAME) where FILENAME has N paused segments
%    will create N individual files with the same extension named
%    FILENAME-1, FILENAME-2, ..., FILENAME-N.
%
%    separatePausedFiles without any input arguments opens a UIgetfile to
%    select the NSx file to separate
%
%    Brett Dowden
%    bdowden@blackrockmicro.com
%    Nick Halper
%    nhalper@blackrockmicro.com
%    Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0: Initial release
%
% 1.1.0.0:
%   - Fixed a broken function that dependent on a non-existant saveNEV
%   script.
%
% 1.2.0.0:
%   - Corrected  NSx_out.MetaTags.Filename     =
%   [NSx.MetaTags.Filename(1:end) '-p' sprintf('%03d', i)
%   NSx.MetaTags.FileExt] which was incorrectly removing the last 4
%   characters of the filename in an attempt to remove the extension.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since openNSx checks input parameters for file validity and for no input
% argument case, just pass varargin to openNSx
NSx = openNSx('read',varargin{:});

% Check to see if there were pauses in the file.  If so, save the sections
% as individual files, if not just report that to the user and end function
if length(NSx.MetaTags.Timestamp) == 1
    disp('No pauses in data file.  No action taken');
else
    for i = 1 : length(NSx.MetaTags.Timestamp)
        % create a copy of NSx, except with only one data segment and the
        % corresponding Timestamp, DataPoints, and NumofPacket fields
        NSx_out.MetaTags              = NSx.MetaTags;
        NSx_out.MetaTags.Timestamp    = NSx.MetaTags.Timestamp(i);
        NSx_out.MetaTags.DataPoints   = NSx.MetaTags.DataPoints(i);
        NSx_out.MetaTags.NumofPackets = NSx.MetaTags.DataPoints(i);
        
        NSx_out.ElectrodesInfo        = NSx.ElectrodesInfo;
        
        NSx_out.Data                  = NSx.Data{i};
        NSx_out.RawData.Headers       = NSx.RawData.Headers;
        %Data headers are only 9 bytes, so only check for the correct,
        %respective 9 bytes
        NSx_out.RawData.DataHeader    = NSx.RawData.DataHeader(1+(9*(i-1)):9+(9*(i-1)));
        
        % Create enumerated file name
        NSx_out.MetaTags.Filename     = [NSx.MetaTags.Filename(1:end) '-p' sprintf('%03d', i) NSx.MetaTags.FileExt];
        disp('Opening the original file...');


% Check for filename existence
newFilename = fullfile(NSx.MetaTags.FilePath, NSx_out.MetaTags.Filename);
if exist(newFilename, 'file') == 2
    overwriteFlag = input('The file already exists. Overwrite? ', 's');
    if ~strcmpi(overwriteFlag, 'y')
        return;
    end
end
        

% Create new file
FIDw = fopen(newFilename, 'w+', 'ieee-le');


% Writing the header information
fwrite(FIDw, NSx_out.RawData.Headers, 'uint8');
fwrite(FIDw, NSx_out.RawData.DataHeader, 'uint8');


% Writing data into file
disp('Writing into the new file...');
fwrite(FIDw, NSx_out.Data, 'int16');
fclose(FIDw);
    end
    clear all;
end