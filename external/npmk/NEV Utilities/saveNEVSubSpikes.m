function saveNEVSubSpikes(channelsToRead, fileFullPath, addedSuffix)

% saveNEVSubSpikes
% 
% Opens saves a new NEV file that only contains chanenls in channelsToRead.
%
% Use saveNEVSubSpikes(channelsToRead)
%
%   channelsToRead: The channel data to be saved into a new NEV.
%                   DEFAULT: This input is required.
%
%   fileFullPath:   The full path to the NEV file that is going to be
%                   used for splitting.
%                   DEFAULT: The user will be prompted to choose a suffix.
%
%   addedSuffix:    The suffix added to the end of splitted file.
%                   DEFAULT: 'ss' is the added suffix.
%
%   Example 1:
%   channelsToRead(4, 'c:\datafolder\datafile.nev', 'tet');
%
%   In the example above, the file c:\datafolder\datafile.nev will be used.
%   The selected NEV file will be saved into a new NEV file that only
%   contains data from channel 4. The new file will have the added suffix
%   'tet', so the new filename will be c:\datafolder\datafile-ssXXX.nev.
%
%   Example 2:
%   channelsToRead([5,8,12];
%
%   In the example above, the user will be prompted to select a NEV file.
%   The selected NEV file will be saved into a new NEV file that only
%   contains data from channels 5, 8, and 12.
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   Version 1.2.2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0: Initial release.
% 
% 1.1.0.0:
%   - Bug fix with saving file names.
%   - Added ability to pass the file name to the function as an argument.
%   - Added ability to define the suffix to the split files (addedExtension).
%
% 1.2.0.0:
%   - Fixed a bug related to the # of input arguments and compatibility
%     with other functions.
%
% 1.2.1.0:
%   - Updated help.
%
% 1.2.2.0:
%   - Fixed a bug where the data was not being saved correctly on Windows
%     machines.
%   - Fixed a bug where tetrodes higher than 10 were overwriting tetrodes 1
%     through 10 over and over again.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Validating the input argument
if ~exist('channelsToRead', 'var')
    disp('channelsToRead is a required input argument.');
    return;
else
    channelsToRead = [channelsToRead 65535:-1:65520];
end

%% Setting the default for addedExtension
if ~exist('addedSuffix', 'var')
    addedSuffix = 'ss';
end

%% Opening the file and reading header
if ~exist('fileFullPath', 'var')
    [dataFilename dataFolder] = getFile('*.nev');
    fileFullPath = [dataFolder dataFilename];
end
FID = fopen(fileFullPath, 'r', 'ieee-le');

%% Calculating the header bytes
BasicHeader = fread(FID, 20, '*uint8');
headerBytes  = double(typecast(BasicHeader(13:16), 'uint32'));
dataPacketByteLength = double(typecast(BasicHeader(17:20), 'uint32'));

%% Calculating the data file length and eeking to the beginning of the file
fseek(FID, 0, 'eof');
endOfDataByte = ftell(FID);
dataByteLength = endOfDataByte - headerBytes;
numberOfPackets = double(dataByteLength) / double(dataPacketByteLength);

%% Reading the header binaries and saving it for future
fseek(FID, 0, 'bof');
headerBinaries = fread(FID, headerBytes, '*uint8');

%% Reading the data binaries
dataBinaries = fread(FID, [dataPacketByteLength numberOfPackets], '*uint8', 0);

%% Finding what PacketIDs have the desired channels
for IDX = 1:size(dataBinaries,2)
    PacketIDs(IDX) = typecast(dataBinaries(5:6, IDX), 'uint16');
end
packetIDIndices = [];
for IDX = 1:length(channelsToRead)
    packetIDsFound = find(PacketIDs == channelsToRead(IDX));
    packetIDIndices(length(packetIDIndices)+1:length(packetIDIndices)+length(packetIDsFound)) = packetIDsFound;
end
packetIDIndices = sort(packetIDIndices);

%% Truncating the data to only contain the desired channels
newDataBinaries = dataBinaries(:, packetIDIndices);

%% Determining the file name
currentFileNames = dir([fileFullPath(1:end-4) '-' addedSuffix '*.nev']);
if ~isempty(currentFileNames)
    fileIDX = str2num(currentFileNames(end).name(end-6:end-4))+1;
else
    fileIDX = 1;
end

%% Saving the new NEV containig the desired channels
FIDw = fopen([fileFullPath(1:end-4) '-' addedSuffix sprintf('%03d', fileIDX) fileFullPath(end-3:end)], 'w+', 'n');
fwrite(FIDw, headerBinaries, 'uint8');
fwrite(FIDw, newDataBinaries, 'uint8');
fclose(FID);
fclose(FIDw);