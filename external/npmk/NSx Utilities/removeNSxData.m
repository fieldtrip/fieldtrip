function removeNSxData(timestamps, filename)
%
%   removeNSxData(timestamps, filename)
%
%   This function removes chunks of data given in a 2-column timestamps
%   variable from a NSx file and then saves it as another NSx. This is
%   useful when removing useless data from the file to decrease the size of
%   the NSx file without modifying the file timestamps. The file size only
%   changes when the NSx file is compressed.
%
%   timestamps: This variable contains beginning and ending timestamps of
%               the data that need to be removed from the NSx file. The
%               format is:
%                          [BegTimeStamp1 EndTimeStamp1
%                           BegTimeStamp2 EndTimeStamp2
%                           BegTimeStamp3 EndTimeStamp3
%                                ...          ...
%                           BegTimeStampN EndTimeStampN]
%
%               The data in the NSx between the pairs of timestamps will be
%               removed and set to 0.
%               (REQUIRED)
%
%   filename:   The name of the NSx file to be used. This input is also
%               optional. In its absense, a dialog will open and will
%               prompt the user to select an NSx file.
%               (OPTIONAL)
%
%   Example:    removeNSxData([100,600; 1000,2000], 'c:\datafile\sampleNSx.ns5')
%
%               In the example above, the timestamps between 100:600 and
%               1000:2000 will be removed from the file sampleNSx.ns5 and
%               the resulting NSx will be saved in sampleNSx_chunked.ns5
%               file.
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%
%   Version 1.0.1.0

% Validating the timestamps variable
if ~exist('timestamps', 'var')
    disp('The variable timestamps is required.');
    return;
end

if size(timestamps, 2) ~= 2
    disp('The timestamps variable should be a Nx2 variable containing beginning timestamps in the first column and ending timestamps in the second column.');
    return;
end

if ~all(timestamps(:,1) < timestamps(:,2))
    disp('The timestamps in the first column (beginning) have to be smaller than the ones in the second column (ending).')
    return;
end
    
% Figuring out the filename, if not passed on to the function
if ~exist('filename', 'var')
	[fname, path] = getFile('*.*', 'Choose an NSx file...');
    filename = [path fname];
end

if exist(filename, 'file') ~= 2
    disp('The file does not exist.');
    return;
end

% Openning and reading the file
FID = fopen(filename, 'r', 'ieee-le');
readData = fread(FID, '*uint8');
FID = fclose(FID);

% Extracting some information from the NSx file
fileTypeID   = char(readData(1:8))';
headerBytes  = typecast(readData(11:14), 'uint32')+9;
samplingFreq = double(typecast(readData(291:294), 'uint32'))/...
               double(typecast(readData(287:290), 'uint32'));
channelCount = typecast(readData(311:314), 'uint32');

% Converting timestamps into uint16 (2 bytes per data point)
timestamps = timestamps*2;

% Removing given timestamps
for idx = 1:size(timestamps, 1)
    readData((timestamps(idx,1)-2)*channelCount+1+headerBytes:timestamps(idx,2)*channelCount+headerBytes) = NaN;
end

% Saving the chunked file
writeFilename = [filename(1:end-4) '_chunked', filename(end-3:end)];
FID = fopen(writeFilename, 'w+', 'ieee-le');
fwrite(FID, readData);
fclose(FID);    