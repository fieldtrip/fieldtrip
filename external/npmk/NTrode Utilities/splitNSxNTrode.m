function splitNSxNTrode

% splitNSxNTrode
% 
% Opens and splits an NSx file based on its NTrode groups. It depends on
% openCCF file.
%
% Use splitNSxNTrode
% 
% This function does not take any inputs.
%
%   Example 1: 
%   splitNSxNTrode;
%
%   In the example above, the user will be prompted to select a CCF file
%   first. The CCF contains the ntrode grouping infromation. Then the user
%   will be prompted to select a NSx file. The script will then split the
%   NSx file into smaller NSx files containing channels in given ntrode
%   groups. For example, if ntrode group one consists of channels 1,3,5,
%   and 12, then using splitNSxNTrode will split the file into a smaller
%   NSx file that contains those channels only. If there are multiple
%   ntrodes then the files will split into multiple smaller files, equal in
%   number of the ntrode groups.
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   Version 1.0.1.0
%

% Validating input parameter
ccf = openCCF;
splitCount = length(ccf.NTrodeInfo.NTrodeID);

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
    positionEOB = ftell(FID);
    positionEOE = typecast(basicHeader(11:14), 'uint32');
    fseek(FID, 0, 'eof');
    positionEOD = ftell(FID);
    
    % Calculating channelCount, data Length
    channelCount = typecast(basicHeader(311:314), 'uint32');
    dataLength = positionEOD - positionEOE - 9;
    
    % Reading extended header and calculating channel IDs
    fseek(FID, positionEOB, 'bof');
    extHeader = fread(FID, [(positionEOE-positionEOB)/channelCount, channelCount], '*uint8');
 
    % Reading the channel IDs. This only wokrs for as long as channel IDs
    % are 8-bit integers. If they change, a typecast of 3:4 is necessasry
    channelID = extHeader(3,1:channelCount);
    
    % Reading the number of packets
    fseek(FID, 28, 'bof');
   
    % Calculating the number of bytes in each segment
    channelBytes = (dataLength)/channelCount;

    % Reading the headers and the data header
    fseek(FID, 0, 'bof');
    fileHeader = fread(FID, positionEOE, 'char');
    dataHeader = fread(FID, 9, 'char');
	fseek(FID, positionEOE+9, 'bof');

    % Reading the data
    fprintf('\nReading the entire data file...\n');
    dataSegment = fread(FID, [channelCount, channelBytes], 'int16');
       
    for idx = 1:splitCount
        % Determining whether tetrode channel is recorded and valid in NSx
        tetrodeChannels = ccf.NTrodeInfo.NTrodeMembers{idx};
        for tetIDX = 1:length(tetrodeChannels)
            validChannel = find(channelID == tetrodeChannels(tetIDX)); 
            if isempty(validChannel)
                fprintf(2,'The tetrode channel %1.0f from tetrode group %d does not exist in the continuous file.\n',tetrodeChannels(tetIDX), idx);
                break;
            end
            validChannels(tetIDX) = validChannel;
        end
        if ~isempty(validChannel)
            % Opening a file for saving
            fprintf('Writing segment %d...\n', idx);
            FIDw = fopen([path fname(1:end-4) '-tet' sprintf('%03d', idx) fname(end-3:end)], 'w+', 'ieee-le');
        
            % Writing the segmented data into file
            basicHeader(end-3) = length(validChannels);
            fwrite(FIDw, basicHeader, 'char');
            fwrite(FIDw, extHeader(:,validChannels), 'char');
            fwrite(FIDw, dataHeader, 'char');
            fwrite(FIDw, dataSegment(validChannels,:), 'int16');
            fclose(FIDw);
        end
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
