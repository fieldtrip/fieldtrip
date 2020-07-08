function splitNSx(splitCount)

% splitNSx
% 
% Opens and splits an NSx file in smaller pieces, timewise.
%
% Use splitNSx(splitCount)
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   splitCount:   Defines the number of splits.
%                 DEFAULT: Splits the file in 2 pieces.
%
%   Example 1: 
%   splitNSx(4);
%
%   In the example above, the user will be prompted to select a file. The
%   loaded file will be split in 4 samller files. For example, if the file
%   is 1 hour long then it will be split into four 15-minute files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version History
%
% 1.0.0.0:
%   - Initial release.
%
% 1.1.0.0:
%   - Fixed a bug related to a case where initial timestamp of the first
%     data segment was not 0. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Validating input parameter
if ~exist('splitCount', 'var')
    splitCount = 2;
end

% Getting the file name
[fname, path] = getFile('*.nev', 'Choose an NSx file...');
if fname == 0
    disp('No file was selected.');
    if nargout
        clear variables;
    end
    return;
end
fext = fname(end-3:end);



% Loading the file
%% Reading Basic Header from file into NSx structure.
FID                       = fopen([path fname], 'r', 'ieee-le');
BasicHeader               = fread(FID, 336, '*uint8');
Trackers.fExtendedHeader  = double(typecast(BasicHeader(13:16), 'uint32'));
Trackers.countPacketBytes = double(typecast(BasicHeader(17:20), 'uint32'));

%% Doing Trackers
fseek(FID, 0, 'eof');
Trackers.fData = ftell(FID);
Trackers.countDataPacket = (Trackers.fData - Trackers.fExtendedHeader)/Trackers.countPacketBytes;


fseek(FID, Trackers.fExtendedHeader, 'bof');
tRawData  = fread(FID, [10 Trackers.countDataPacket], '10*uint8=>uint8', Trackers.countPacketBytes - 10);
Timestamp = tRawData(1:4,:);
Timestamp = typecast(Timestamp(:), 'uint32').';

splitPacketStarts = find(diff(Timestamp)<0);
splitPacketBytes = Trackers.countPacketBytes * splitPacketStarts;

% Reading headers and seeking to beginning of data
fseek(FID, 0, 'bof');
fileHeader = fread(FID, Trackers.fExtendedHeader, '*uint8');

for idx = 1:length(splitPacketStarts)
    % Opening a file for saving
    FIDw = fopen([path fname(1:end-4) '-s' sprintf('%03d', idx) fname(end-3:end)], 'w+', 'ieee-le');
    fprintf('\nReading segment %d... ', idx);
    % Reading the segment
    dataSegment = fread(FID, splitPacketBytes(idx), 'char');
    fprintf('Writing segment %d... ', idx);
    % Writing the segmented data into file
    fwrite(FIDw, fileHeader, 'char');
    fwrite(FIDw, dataSegment, 'char');
    % Clearing variables and closing file
    clear dataSegment;
    fclose(FIDw);
end