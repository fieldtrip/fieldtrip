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
if ~ismac
    [fname, path] = getFile('*.ns*', 'Choose an NSx file...');
else
    [fname, path] = getFile('*.*', 'Choose an NSx file...');
end
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
NSx.MetaTags.Filename     = fname;
NSx.MetaTags.FilePath     = path(1:end-1);
NSx.MetaTags.FileExt      = fext;
NSx.MetaTags.FileTypeID   = fread(FID, [1,8]   , '*char');
disp(['Splitting the NSx file in ' num2str(splitCount) ' pieces...']);
if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALSG')
    disp('File type 2.1 is not yet implemented.');
    %NOT IMPLEMENTED YET
%     fseek(FID, 0, 'bof');
%     header = fread(FID, 314,'*uint8');
%     positionEOH = ftell(FID);
%     fseek(FID, 0, 'eof');
%     positionEOD = ftell(FID);
%     dataLength = positionEOD - positionEOH;
%     fseek(FID, 28, 'bof');
%     channelCount = fread(FID, 1       , 'uint32=>double');
elseif strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD')
    % Calculating different points in the file
    fseek(FID, 0, 'bof');
    basicHeader = fread(FID, 314, '*uint8');
    positionEOE = typecast(basicHeader(11:14), 'uint32');
    fseek(FID, 0, 'eof');
    positionEOD = ftell(FID);
    % Calculating channelCount, data Length
    channelCount = typecast(basicHeader(311:314), 'uint32');
    dataLength = positionEOD - positionEOE - 9;
    % Reading the number of packets
    fseek(FID, 28, 'bof');
    numOfPackets = (dataLength)/(2*channelCount);
    % Calculating the number of bytes in each segment
    segmentBytes = floor(numOfPackets/splitCount)*(2*channelCount);
    % Reading the headers and the data header
    fseek(FID, 0, 'bof');
    fileHeader = fread(FID, positionEOE, 'char');
    dataHeader = fread(FID, 9, 'char');
	fseek(FID, positionEOE+9, 'bof');
    for idx = 1:splitCount
        % Opening a file for saving
        FIDw = fopen([path fname(1:end-4) '-s' sprintf('%03d', idx) fname(end-3:end)], 'w+', 'ieee-le');
        fprintf('\nReading segment %d... ', idx);
        % Reading the segment
        dataSegment = fread(FID, segmentBytes, 'char');
        fprintf('Writing segment %d... ', idx);
        % Writing the segmented data into file
        fwrite(FIDw, fileHeader, 'char');
        % Set the timestamp of the segments 2+ to 0 so there's no
        % introduced shift by openNSx.
        if idx > 1
            dataHeader(2:5) = 0;
        end
        fwrite(FIDw, dataHeader, 'char');
        fwrite(FIDw, dataSegment, 'char');
        % Clearing variables and closing file
        clear dataSegment;
        fclose(FIDw);
    end
    fprintf('\n');
else
    % Display error if non-compatible file is trying to open.
    disp('This version of splitNSx can only split File Specs 2.2 and 2.3');
    disp(['The selected file spec is ' NSx.MetaTags.FileSpec '.']);
    fclose(FID);
    clear variables;
    return;
end
