function splitNSxPauses(fileName)

% splitNSxPauses
% 
% Opens and splits an NSx file in smaller pieces, timewise.
%
% Use splitNSx(fileName)
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   fileName:   File name of the file that needs to be split.
%               DEFAULT: The user will be prompted to select a file.
%
%   Example 1: 
%   splitNSx('C:\Datafolder\mydata.ns5');
%
%   In the example above, the file C:\Datafolder\mydata.ns5 will be opened.
%   The loaded file will be split in samller files representing its paused 
%   sub-segments.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version History
%
% 1.0.0.0: August 31, 2016
%   - Initial release.
%   - Successor to separateNSxPaused running much more memory efficient.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Getting the file name
if ~exist('fileName', 'var')
    fileName = '';
end

if ~(exist(fileName, 'file') == 2)
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
else
    [path, fname, fext] = fileparts(fileName);
    if ismac path = [path '/']; else path = [path '\']; end
    fname = [fname fext];
end

%% Getting header information
NSx = openNSx('noread', [path fname ]);
    
% Loading the file
%% Reading Basic Header from file into NSx structure.
FID                       = fopen([path fname], 'r', 'ieee-le');
NSx.MetaTags.Filename     = fname;
NSx.MetaTags.FilePath     = path(1:end-1);
NSx.MetaTags.FileExt      = fext;
NSx.MetaTags.FileTypeID   = fread(FID, [1,8]   , '*char');
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
    % Calculating the number of splits
    splitCount = length(NSx.MetaTags.Timestamp);
    % Calculating the number of bytes in each segment
    segmentBytes = NSx.MetaTags.DataPoints * 2 * double(channelCount);
    % Reading the headers and the data header
    fseek(FID, 0, 'bof');
    fileHeader = fread(FID, positionEOE, 'char');
    dataHeader = fread(FID, 9, 'char');
	fseek(FID, positionEOE, 'bof');
    disp(['Splitting the NSx file in ' num2str(splitCount) ' pieces...']);
    for idx = 1:splitCount
        % Opening a file for saving
        FIDw = fopen([path fname(1:end-4) '-s' sprintf('%03d', idx) fname(end-3:end)], 'w+', 'ieee-le');
        fprintf('\nReading segment %d... ', idx);
        % Reading the segment
        fseek(FID, 9, 'cof'); % Skipping the data header
        dataSegment = fread(FID, segmentBytes(idx), 'char');
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
