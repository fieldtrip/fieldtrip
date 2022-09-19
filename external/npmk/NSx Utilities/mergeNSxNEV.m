%
%   mergeNSxNEV(filename1, filename2)
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
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%
%   Version 1.2.2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0: Initial release
%
% 1.2.2.0: August 3, 2016
%   - Fixed a bug that resulted in a crash if one of two NEV files weren't
%     available.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mergeNSxNEV(varargin)

if nargin == 3
    if strcmpi(varargin{1}, 'nowarning')
        warningMessage = varargin{1};
        filename1 = varargin{2};
        filename2 = varargin{3};
    elseif strcmpi(varargin{3}, 'nowarning')
        filename1 = varargin{1};
        filename2 = varargin{2};
        warningMessage = varargin{3};
    else
        disp('Invalid arguments. See help by typing ''help mergeNSxNEV''');
        return;
    end
elseif nargin == 2
    filename1 = varargin{1};
    filename2 = varargin{2};
    warningMessage = 'warning';
elseif nargin == 1
    warningMessage = varargin{1};
elseif nargin == 0
    warningMessage = 'warning';
else
    disp('Invalid number of arguments. See help by typing ''help mergeNSxNEV''');
    return;
end

%% Warning
if ~strcmpi(warningMessage, 'nowarning')
    fprintf(2, '\nCaution! This scripts assumes that the two segments of NSx ');
    fprintf(2, '\nand NEV files are from the same recording session. The settings ');
    fprintf(2, '\nmust be identical. Always verify the merge to assure it was done');
    fprintf(2, '\ncorrectly. To suppress this message, run mergeNSxNEV with command');
    fprintf(2, '\n''nowarning''. Ex: mergeNSxNEV(''nowarning'').\n\n');
end

%% Getting the file names
if exist('filename1', 'var') && exist('filename2', 'var')
    if exist(filename1, 'file') && exist(filename2, 'file')
        NSx1FullName = filename1;
        NEV1FullName = [NSx1FullName(1:end-3) 'nev'];
        NSx2FullName = filename2;
        NEV2FullName = [NSx2FullName(1:end-3) 'nev'];
    else
        disp('One or both files do not exist.');
        return;
    end
else
    disp('Load the first NSx file.');
    if ~ismac
        [NSx1dataFilename, NSx1dataFolder] = getFile('*.ns*', 'Choose the first NSx file...');
    else
        [NSx1dataFilename, NSx1dataFolder] = getFile('*.*', 'Choose the first NSx file...');
    end
    if NSx1dataFilename == 0
        disp('No file was selected.');
        return;
    end
    disp('Load the second NSx file.');
    if ~ismac
        [NSx2dataFilename, NSx2dataFolder] = getFile('*.ns*', 'Choose the second NSx file...');
    else
        [NSx2dataFilename, NSx2dataFolder] = getFile('*.*', 'Choose the second NSx file...');
    end
    if NSx2dataFilename == 0
        disp('No file was selected.');
        return;
    end
    NSx1FullName = strcat(NSx1dataFolder, NSx1dataFilename);
    NSx2FullName = strcat(NSx2dataFolder, NSx2dataFilename);
	NEV1FullName = [NSx1FullName(1:end-3) 'nev'];
	NEV2FullName = [NSx2FullName(1:end-3) 'nev'];
end

%% Reading NSx1
NSx1 = openNSx(NSx1FullName);
FID1 = fopen(NSx1FullName, 'r', 'ieee-le');
fread(FID1, 10, '*int8');
NSx1HeaderBytes = typecast(NSx1.RawData.Headers(11:14), 'uint32');
fread(FID1, NSx1HeaderBytes+5, '*int8');
NSx1DataLength = NSx1.MetaTags.DataPoints;
fseek(FID1, 0, 'bof');
NSx1Data = fread(FID1, '*int8');
fclose(FID1);

NSx1LengthIn30k = NSx1DataLength * (NSx1.MetaTags.TimeRes / NSx1.MetaTags.SamplingFreq);


%% Reading NSx2
FID2 = fopen(NSx2FullName, 'r', 'ieee-le');
fread(FID2, 10, '*int8');
NSx2HeaderBytes = fread(FID2, 1, '*int32');
fseek(FID2, 0, 'bof');
fread(FID2, NSx2HeaderBytes+5, '*int8');
NSx2DataLength = fread(FID2, 1, '*uint32');
NSx2Data = fread(FID2, '*int8');
fclose(FID2);

%% Combining the two files and adjusting data length
NSx1Data(NSx1HeaderBytes + 5 + 1:NSx1HeaderBytes + 9) = typecast((NSx1DataLength + NSx2DataLength), 'int8');
NSx1Data = [NSx1Data; NSx2Data];


%% Writing data back to combined file

% Writing header into the file
newFilename = [NSx1FullName(1:end-4) '-combined.' NSx1FullName(end-2:end)];
if exist(newFilename, 'file') == 2
    overwriteFlag = input('The NSx file already exists. Overwrite? ', 's');
    if ~strcmpi(overwriteFlag, 'y')
        clear all;
        return;
    end
end

% Writing combined NSx
FIDw = fopen(newFilename, 'w+', 'ieee-le');
disp('Writing NSx data back into the combined file...');
fwrite(FIDw, NSx1Data, 'int8');
fclose(FIDw);

%% Opening the NEV files

if ~exist('NEV1FullName', 'var')
    disp('Load the first NEV file.');
    [NEV1dataFilename NEV1dataFolder] = getFile('*.*');
    disp('Load the second NEV file.');
    [NEV2dataFilename NEV2dataFolder] = getFile('*.*');

    NEV1FullName = strcat(NEV1dataFolder, NEV1dataFilename);
    NEV2FullName = strcat(NEV2dataFolder, NEV2dataFilename);
end

%% Reading NEV1
FID1 = fopen(NEV1FullName, 'r', 'ieee-le');
FID2 = fopen(NEV2FullName, 'r', 'ieee-le');
okNEV = (FID1 ~= -1) && (FID2 ~= -1);
if okNEV
    NEV1Data = fread(FID1, '*int8');
    fclose(FID1);

    %% Reading NEV2
    fseek(FID2, 12, 'bof');
    NEV2HeaderBytes = fread(FID2, 1, '*int32');
    NEV2BytesinDataPackets = fread(FID2, 1, '*int32');
    fseek(FID2, NEV2HeaderBytes, 'bof');
    NEV2Data = fread(FID2, '*int8');
    fclose(FID2);

    %% Adding the timestamp of the first file to the second file
    NEV2Data = typecast(NEV2Data, 'uint32');
    NEV2Data = reshape(NEV2Data, NEV2BytesinDataPackets/4, length(NEV2Data)/(NEV2BytesinDataPackets/4));
    NEV2Data(1,:) = NEV2Data(1,:) + NSx1LengthIn30k;
    NEV2Data = reshape(NEV2Data, size(NEV2Data, 1) * size(NEV2Data, 2), 1);
    NEV2Data = typecast(NEV2Data, 'int8');

    NEV3Data = [NEV1Data; NEV2Data];

    %% Writing data back to combined NEV file

    % Writing header into the file
    newFilename = [NEV1FullName(1:end-4) '-combined.' NEV1FullName(end-2:end)];
    if exist(newFilename, 'file') == 2
        overwriteFlag = input('The NEV file already exists. Overwrite? ', 's');
        if ~strcmpi(overwriteFlag, 'y')
            clear all;
            return;
        end
    end

    % Writing combined NEV
    FIDCombined = fopen(newFilename, 'w+', 'ieee-le');
    disp('Writing NEV data back into the combined file...');
    fwrite(FIDCombined, NEV3Data, 'int8');
    fclose(FIDCombined);
    clear all;
else
    disp('One or more NEV files were not found. Skipping NEV merging.');
end
