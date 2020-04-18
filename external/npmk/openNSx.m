function varargout = openNSx(varargin)

% openNSx
% 
% Opens and reads an NSx file then returns all file information in a NSx
% structure. Works with File Spec 2.1, 2.2, 2.3, and 3.0.
% Use OUTPUT = openNSx(fname, 'read', 'report', 'e:xx:xx', 'c:xx:xx', 't:xx:xx', 'mode', 'precision', 'skipfactor').
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   fname:        Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file. 
%                 DEFAULT: Will open Open File UI.
%
%   'read':       Will read the data in addition to the header information
%                 if user passes this argument.
%                 DEFAULT: will read the entire file.
%
%   'report':     Will show a summary report if user passes this argument.
%                 DEFAULT: will not show report.
%
%   'e:XX:YY':    User can specify which electrodes need to be read. The
%                 number of electrodes can be greater than or equal to 1
%                 and less than or equal to 128. The electrodes can be
%                 selected either by specifying a range (e.g. 20:45) or by
%                 indicating individual electrodes (e.g. 3,6,7,90) or both.
%                 Note that, when individual channels are to be read, all
%                 channels in between will also be read. The prorgam will
%                 then remove the unwanted channels. This may result in a
%                 large memory footprint. If memory issues arrise please
%                 consider placing openNSx in a for loop and reading
%                 individual channels.
%                 This field needs to be preceded by the prefix 'e:'. See
%                 example for more details. If this option is selected the
%                 user will be promped for a CMP mapfile (see: KTUEAMapFile)
%                 provided by Blackrock Microsystems. This feature required
%                 KTUEAMapFile to be present in path.
%                 DEFAULT: will read all existing electrodes.
%
%   'c:XX:YY':    User can specify which channels need to be read. The
%                 number of channels can be greater than or equal to 1
%                 and less than or equal to 255. The channels can be
%                 selected either by specifying a range (e.g. 20:45) or by
%                 indicating individual channels (e.g. 3,6,7,90) or both.
%                 Note that, when individual channels are to be read, all
%                 channels in between will also be read. The prorgam will
%                 then remove the unwanted channels. This may result in a
%                 large memory footprint. If memory issues arrise please
%                 consider placing openNSx in a for loop and reading
%                 individual channels.
%                 This field needs to be preceded by the prefix 'c:'. See
%                 example for more details.
%                 DEFAULT: will read all existing analog channels.
%
%   't:XX:YY':    User can specify the beginning and end of the data
%                 segment to be read. If the start time is greater than the
%                 length of data the program will exit with an errorNS
%                 message. If the end time is greater than the length of
%                 data the end packet will be selected for end of data. The
%                 user can specify the start and end values by comma 
%                 (e.g. [20,50]) or by a colon (e.g. [20:50]). To use this
%                 argument the user must specify the [electrodes] or the
%                 interval will be used for [electrodes] automatically.
%                 This field needs to be preceded by the prefix 't:'. 
%                 Note that if 'mode' is 'sample' the start duration cannot
%                 be less than 1. The duration is inclusive.
%                 See example for more details.
%                 DEFAULT: will read the entire file.
%
%   'mode':       The user can specify the mode of duration in [duration],
%                 such as 'sec', 'min', 'hour', or 'sample'. If 'sec' is
%                 specified the numbers in [duration] will correspond to
%                 the number of seconds. The same is true for 'min', 'hour'
%                 and 'sample'.
%                 DEFAULT: reads 'sample'.
%
%   'uV':         Will read the spike waveforms in unit of uV instead of
%                 raw values. Note that this conversion may lead to loss of
%                 information (e.g. 15/4 = 4) since the waveforms type will
%                 stay in int16. It's recommended to read raw spike
%                 waveforms and then perform the conversion at a later
%                 time.
%                 DEFAULT: will read waveform information in raw.
%
%   'precision':  This will specify the precision for NSx file. If set to
%                 'double' the NSx data will be read as 'double' and if set
%                 to 'short', the NSx data will be read as 'int16' data
%                 type. While reading the file as 'short' may have a much
%                 smaller memory footprint and a faster read time, some 
%                 post data analysis such as multiplying the signal by a 
%                 factor that will make the data larger than (-32,768 to 
%                 32,767 -- refer to MATLAB documentation for more 
%                 information) may result in unexpected behavior. 
%                 Always use caution when using short. If you are not sure
%                 of what to use then do not specify this option.
%                 DEFAULT: will read data in 'int16'.
%
%   'skipfactor': This option will allow the user to read a decimated
%                 version of the data. The skipfactor will determine how
%                 many samples to skip. For example, if skipfactor is 2
%                 then every other sample is read. If skipfactor is 5 then
%                 every fifth sample is read. This is useful to briefly
%                 looking at the data in a large datafile when reading the
%                 entire dataset would overflow the memory.
%                 DEFAULT: is set to 1, so every sample will be read.
%
%   'ver':        If this argument is passed to the function it will return
%                 the version number of the function without reading any
%                 data files.
%
%   OUTPUT:       Contains the NSx structure.
%
%   Example 1: 
%   openNSx('report','read','c:\data\sample.ns5', 'e:15:30', 't:3:10','min', 'p:short', 's:5');
%
%   or equivalently
%   openNSx('report','read','c:\data\sample.ns5', 'electrodes', 15:30, 'duration', 3:10, 'min', 'precision', 'short', 'skipfactor', 5);
%
%   In the example above, the file c:\data\sample.ns5 will be used. A
%   report of the file contents will be shown. The data will be read from
%   electrodes 15 through 50 in the 3-10 minute time interval. A decimated 
%   version of the datafile will be read, where only every 5th sample is read.
%   If any of the arguments above are omitted the default values will be used.
%
%   Example 2:
%   openNSx('read','c:15:30');
%
%   In the example above, the file user will be prompted for the file. The
%   file will be read using 'int16' precision as default. All time points 
%   of Only channels 15 through 30 will be read. If any of the arguments 
%   above are omitted the default values will be used.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Version History
%
% 5.1.8.2:
%   - Fixed the way DayOfWeek is read in MetaTags.
%
% 5.1.9.0:
%   - Fixed a bug where with skipFactor being read correctly as a num.
%
% 5.1.10.0:
%   - Updated feature to save data headers for a paused file. It is a
%     dependent feature for seperatePausedNSx.
%
% 5.1.11.0:
%   - Fixed an issue where 1 sample would not be read when using the
%     t:xx:xx argument and 'sample'.
%   - Fixed an error when 'duration' was used to load specific data length.
%
% 5.1.12.0:
%   - Better error handling if a file is not provided and an output
%     variable was requested by the calling function.
%
% 5.2.0.0: June 12, 2014
%   - It removes the extra ElectrodesInfo entried for channels not
%     read if 'c:XX:XX' or 'e:XX:XX' are used.
%   - It reports variable ChannelCount under MetaTags correctly.
%   - It automatically compensate for any NSx file with non-0 beginnings
%     and adds 0s for to the begining of the file to properly align the
%     timestamps.
%
% 5.2.1.0: June 12, 2014
%   - Fixed a small bug where extra 0s were tacked on to the beginning of
%     paused file segments.
%   - Updated the version.
%
% 5.2.2.0: June 13, 2014
%   - Fixed bug for when 'noread' was used on a paused file.
%
% 6.0.1.0: December 2, 2014
%   - Fixed a bug related to file format 2.1 not being read correctly.
%   - Corrected the way Filename, FileExt, and FilePath was being
%     processed.
%   - File dialogue now only shows NSx files on non Windows-based
%     computers.
%   - Added 512 synchronized reading capability.
%   - Now on non-Windows computers only NSx files are shown in the file
%     dialogue.
%   - Fixed the date in NSx.MetaTags.DateTime.
%
% 6.1.0.0: March, 15 2015
%   - Added the ability to read from networked drives in Windows.
%   - Fixed the DateTime variable in MetaTags.
%   - Fixed the date in NSx.MetaTags.DateTime (again).
%   - Fixed a bug related to starting and stopping packets when a specific
%     time is passed to the function.
%   - Fixed a bug where 512+ ch rules were being applied to smaller channel
%     count configuration.
%
% 6.1.1.0: June 15, 2015
%   - Bug fixes related to timestamps when the recording didn't start at
%     proctime 0.
%
% 6.2.0.0: October 1, 2015
%   - Fixed a bug related to reading the correct length of time when a skip
%     factor was used.
%   - Bug fixes related to information that separatePausedNSx depends on.
%   - Added 'uV' as an option to read the data in the unit of uV.
%
% 6.2.1.0: April 16, 2016
%   - Fixed a bug related to converting the unit to uV in case of having
%     multiple data segments (paused file).
%
% 6.2.2.0: July 6, 2016
%   - Fixed another bug related to converting the unit to uV.
%
% 6.3.0.0: August 3, 2016
%   - Added support for loading a segment of paused files.
%
% 6.3.1.0: August 31, 2016
%   - Fixed a bug when reading a non-o start across a paused segment.
%
% 6.4.0.0: December 1, 2016
%   - Fixed a serious bug related to loading paused files.
%   - Fixed a bug where an empty data segment resulted in a cell structure.
%
% 6.4.1.0: June 15, 2017
%   - It is no longer necessary to provide the full path for loading a
%     file.
%
% 6.4.2.0: September 1, 2017
%   - Fixed a bug related to reading data from sample that is not 1 and
%     timestamp that used to get reset to 0.
%
% 6.4.3.0: September 13, 2017
%   - Removed a redundant block of code that was accidentally placed in the
%     script twice.
%   - Checks to see if there's a newer version of NPMK is available.
%
% 6.4.3.1: January 24, 2020
%   - Changed file opening access from r+ to r.
%
% 7.0.0.0: January 27, 2020
%   - Added support for 64-bit timestamps in NEV and NSx.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defining the NSx data structure and sub-branches.
NSx          = struct('MetaTags',[],'Data',[], 'RawData', []);
NSx.MetaTags = struct('FileTypeID',[],'SamplingLabel',[],'ChannelCount',[],'SamplingFreq',[], 'TimeRes', [], ...
                      'ChannelID',[],'DateTime',[],'DateTimeRaw',[], 'Comment', [], 'FileSpec', [], ...
                      'Timestamp', [], 'DataPoints', [], 'DataDurationSec', [], 'openNSxver', [], 'Filename', [], 'FilePath', [], ...
                      'FileExt', []);

                                    
NSx.MetaTags.openNSxver = '7.0.0.0';
                  
%% Check for the latest version fo NPMK
NPMKverChecker

% Defining constants
ExtHeaderLength = 66;
elecReading     = 0;
maxNSPChannels  = 128;
NSx.RawData.PausedFile = 0;
syncShift = 0;

%% Validating the input arguments. Exit with error message if error occurs.
next = '';
for i=1:length(varargin)
    inputArgument = varargin{i};
    if strcmpi(inputArgument, 'ver')
        varargout{1} = NSx.MetaTags.openNSxver;
        return;
    elseif strcmpi(inputArgument, 'channels')
        next = 'channels';
    elseif strcmpi(inputArgument, 'skipfactor')
        next = 'skipfactor';
    elseif strcmpi(inputArgument, 'electrodes')
        next = 'electrodes';
    elseif strcmpi(inputArgument, 'duration')
        next = 'duration';
    elseif strcmpi(inputArgument, 'precision')
        next = 'precision';
    elseif strcmpi(inputArgument, 'report')
        Report = inputArgument;
    elseif strcmpi(inputArgument, 'noread')
        ReadData = inputArgument;
    elseif strcmpi(inputArgument, 'nomultinsp')
        multinsp = 'no';
    elseif strcmpi(inputArgument, 'uV')
        waveformUnits = 'uV';
    elseif strcmpi(inputArgument, 'read')
        ReadData = inputArgument;
    elseif (strncmp(inputArgument, 't:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'duration')
        if strncmp(inputArgument, 't:', 2)
            inputArgument(1:2) = [];
            inputArgument = str2num(inputArgument);
        end
        modifiedTime = 1;
        StartPacket = inputArgument(1);
        EndPacket = inputArgument(end);
        next = '';
    elseif (strncmp(inputArgument, 'e:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'electrodes')
        if exist('KTUEAMapFile', 'file') == 2
            Mapfile = KTUEAMapFile;
            Elec = str2num(inputArgument(3:end)); %#ok<ST2NM>
            if min(Elec)<1 || max(Elec)>128
                disp('The electrode number cannot be less than 1 or greater than 128.');
                if nargout; varargout{1} = -1; end
                return;
            end
            for chanIDX = 1:length(Elec)
                userRequestedChannels(chanIDX) = Mapfile.Electrode2Channel(Elec(chanIDX));
            end
            elecReading = 1;
        else
            disp('To read data by ''electrodes'' the function KTUEAMapFile needs to be in path.');
            clear variables;
            if nargout; varargout{1} = -1; end
            return;
        end
        next = '';
    elseif (strncmp(inputArgument, 's:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'skipFactor')
        if strncmp(inputArgument, 's:', 2)
            skipFactor = str2num(inputArgument(3:end)); %#ok<ST2NM>
        else
            if ischar(inputArgument)
                skipFactor = str2num(inputArgument);
            else
                skipFactor = inputArgument;
            end
        end
        next = '';
    elseif (strncmp(inputArgument, 'c:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'channels')
        if strncmp(inputArgument, 'c:', 2)
            userRequestedChanRow = str2num(inputArgument(3:end)); %#ok<ST2NM>
        else
            userRequestedChanRow = inputArgument;
        end
        next = '';
    elseif (strncmp(varargin{i}, 'p:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'precision')
        if strncmp(varargin{i}, 'p:', 2)
            precisionTypeRaw = varargin{i}(3:end);
        else
            precisionTypeRaw = varargin{i};
        end
        switch precisionTypeRaw
			case 'int16'
				precisionType = '*int16=>int16';
            case 'short'
                precisionType = '*short=>short';
            case 'double'
                precisionType = '*int16';
            otherwise
                disp('Read type is not valid. Refer to ''help'' for more information.');
                if nargout; varargout{1} = -1; end
                return;
        end
        clear precisionTypeRaw;
        next = '';
    elseif strfind(' hour min sec sample ', [' ' inputArgument ' ']) ~= 0
        TimeScale = inputArgument;
    else
        temp = inputArgument;
        if length(temp)>3 && ...
                (strcmpi(temp(3),'\') || ...
                 strcmpi(temp(1),'/') || ...
                 strcmpi(temp(2),'/') || ...
                 strcmpi(temp(1:2), '\\') || ...
                 strcmpi(temp(end-3), '.'))
            fname = inputArgument;
            if exist(fname, 'file') ~= 2
                disp('The file does not exist.');
                if nargout; 
                    varargout{1} = -1; 
                end
                return;
            end
        else
            disp(['Invalid argument ''' inputArgument ''' .']);
            if nargout; varargout{1} = -1; end
            return;
        end
    end
end
clear next;

%% Popup the Open File UI. Also, process the file name, path, and extension
%  for later use, and validate the entry.
if ~exist('fname', 'var')
    [fname, path] = getFile('*.ns1;*.ns2;*.ns3;*.ns4;*.ns5;*.ns6;*.ns6m', 'Choose an NSx file...');
    if fname == 0
        disp('No file was selected.');
        if nargout; varargout{1} = -1; end
        return;
    end
    [~, ~, fext] = fileparts(fname);
else
    if isempty(fileparts(fname))
        fname = which(fname);
    end
    [path,fname, fext] = fileparts(fname);
    fname = [fname fext];
    path  = [path '/'];
end
if fname==0
    if nargout; varargout{1} = -1; end
    return; 
end

%% Loading .x files for multiNSP configuration
if strcmpi(fext(2:4), 'ns6') && length(fext) == 5
    path(1) = fname(end);
    fname(end) = [];
end

tic;

%% Give all input arguments a default value. All input argumens are
%  optional.
if ~exist('Report', 'var');        Report = 'noreport'; end
if ~exist('ReadData', 'var');      ReadData = 'read'; end
if ~exist('StartPacket', 'var');   StartPacket = 1; end
if ~exist('TimeScale', 'var');     TimeScale = 'sample'; end
if ~exist('precisionType', 'var'); precisionType = '*short=>short'; end
if ~exist('skipFactor', 'var');    skipFactor = 1; end
if ~exist('modifiedTime', 'var');  modifiedTime = 0; end
if ~exist('multinsp', 'var');      multinsp = 'yes'; end
if ~exist('waveformUnits', 'var'); waveformUnits = 'raw'; end

% Check to see if 512 setup and calculate offset
if strcmpi(multinsp, 'yes')
    fiveTwelveFlag = regexp(fname, '-i[0123]-');
    if ~isempty(fiveTwelveFlag)
        syncShift = multiNSPSync(fullfile(path, fname));
    else
        multinsp = 'no';
    end
end

if strcmpi(ReadData, 'noread')
%    disp('NOTE: Reading the header information only. To read the data use with parameter ''read'': openNSx(''read'')');
end

if strcmp(Report, 'report')
    disp(['openNSx ' NSx.MetaTags.openNSxver]);
end

%% Reading Basic Header from file into NSx structure.
FID = fopen([path fname], 'r', 'ieee-le');	

fileFullPath = fullfile(path, fname);
[NSx.MetaTags.FilePath, NSx.MetaTags.Filename, NSx.MetaTags.FileExt] = fileparts(fileFullPath);

NSx.MetaTags.FileTypeID   = fread(FID, [1,8]   , '*char');
if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALSG')
	NSx.MetaTags.FileSpec      = '2.1';
    NSx.MetaTags.SamplingLabel = fread(FID, [1,16]  , '*char');
    NSx.MetaTags.TimeRes       = 30000;
    NSx.MetaTags.SamplingFreq  = NSx.MetaTags.TimeRes / fread(FID, 1 , 'uint32=>double');
    ChannelCount               = double(fread(FID, 1       , 'uint32=>double'));
    NSx.MetaTags.ChannelCount  = ChannelCount;
    NSx.MetaTags.ChannelID     = fread(FID, [ChannelCount 1], '*uint32');
    try
    	t                          = dir(fileFullPath);
    	NSx.MetaTags.DateTime      = t.date;
    end
elseif or(strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD'), strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP'))
    BasicHeader                = fread(FID, 306, '*uint8');
    NSx.MetaTags.FileSpec      = [num2str(double(BasicHeader(1))) '.' num2str(double(BasicHeader(2)))];
    HeaderBytes                = double(typecast(BasicHeader(3:6), 'uint32'));
    NSx.MetaTags.SamplingLabel = char(BasicHeader(7:22))';
    NSx.MetaTags.Comment       = char(BasicHeader(23:278))';
    NSx.MetaTags.TimeRes       = double(typecast(BasicHeader(283:286), 'uint32'));
    NSx.MetaTags.SamplingFreq  = NSx.MetaTags.TimeRes / double(typecast(BasicHeader(279:282), 'uint32'));
    t                          = double(typecast(BasicHeader(287:302), 'uint16'));
    ChannelCount               = double(typecast(BasicHeader(303:306), 'uint32'));
    NSx.MetaTags.ChannelCount  = ChannelCount;
    readSize                   = double(ChannelCount * ExtHeaderLength);
    ExtendedHeader             = fread(FID, readSize, '*uint8');
    if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD')
    	timeStampBytes = 4;
    elseif strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP')
        timeStampBytes = 8;
    end
    
    %% Removing extra garbage characters from the Comment field.
    NSx.MetaTags.Comment(find(NSx.MetaTags.Comment==0,1):end) = 0;    
    
	%% Populating extended header information
	for headerIDX = 1:ChannelCount
		offset = double((headerIDX-1)*ExtHeaderLength);
		NSx.ElectrodesInfo(headerIDX).Type = char(ExtendedHeader((1:2)+offset))';
		if (~strcmpi(NSx.ElectrodesInfo(headerIDX).Type, 'CC'))
			disp('extended header not supported');
			fclose(FID);
			if nargout; varargout{1} = -1; end
			return;			
		end
		NSx.ElectrodesInfo(headerIDX).ElectrodeID = typecast(ExtendedHeader((3:4)+offset), 'uint16');
		NSx.ElectrodesInfo(headerIDX).Label = char(ExtendedHeader((5:20)+offset))';
		NSx.ElectrodesInfo(headerIDX).ConnectorBank = char(ExtendedHeader(21+offset) + ('A' - 1));
		NSx.ElectrodesInfo(headerIDX).ConnectorPin   = ExtendedHeader(22+offset);
		NSx.ElectrodesInfo(headerIDX).MinDigiValue   = typecast(ExtendedHeader((23:24)+offset), 'int16');
		NSx.ElectrodesInfo(headerIDX).MaxDigiValue   = typecast(ExtendedHeader((25:26)+offset), 'int16');
		NSx.ElectrodesInfo(headerIDX).MinAnalogValue = typecast(ExtendedHeader((27:28)+offset), 'int16');
		NSx.ElectrodesInfo(headerIDX).MaxAnalogValue = typecast(ExtendedHeader((29:30)+offset), 'int16');
		NSx.ElectrodesInfo(headerIDX).AnalogUnits    = char(ExtendedHeader((31:46)+offset))';
		NSx.ElectrodesInfo(headerIDX).HighFreqCorner = typecast(ExtendedHeader((47:50)+offset), 'uint32');
		NSx.ElectrodesInfo(headerIDX).HighFreqOrder  = typecast(ExtendedHeader((51:54)+offset), 'uint32');
		NSx.ElectrodesInfo(headerIDX).HighFilterType = typecast(ExtendedHeader((55:56)+offset), 'uint16');
		NSx.ElectrodesInfo(headerIDX).LowFreqCorner  = typecast(ExtendedHeader((57:60)+offset), 'uint32');
		NSx.ElectrodesInfo(headerIDX).LowFreqOrder   = typecast(ExtendedHeader((61:64)+offset), 'uint32');
		NSx.ElectrodesInfo(headerIDX).LowFilterType  = typecast(ExtendedHeader((65:66)+offset), 'uint16');
	end
	clear ExtendedHeader;
	%% Parsing and validating FileSpec and DateTime variables
	NSx.MetaTags.DateTimeRaw = t.';
	NSx.MetaTags.DateTime = datestr(datenum(t(1), t(2), t(4), t(5), t(6), t(7)));
	clear t;
else
    disp('This version of openNSx can only read File Specs 2.1, 2.2, 2.3, and 3.0.');
    disp(['The selected file spec is ' NSx.MetaTags.FileSpec '.']);
    fclose(FID);
    clear variables;
	if nargout; varargout{1} = -1; end
    return;
end

% Determining the length of file and storing the value of fEOF
f.EOexH = double(ftell(FID));
fseek(FID, 0, 'eof');
f.EOF = double(ftell(FID));

% Read Raw Header for saveNSx
fseek(FID, 0, 'bof');
NSx.RawData.Headers = fread(FID, f.EOexH, '*uint8');
% if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD')
NSx.RawData.DataHeader = fread(FID, timeStampBytes+5, '*uint8');
% end
fseek(FID, f.EOexH, 'bof');

%% Reading all data headers and calculating all the file pointers for data
% and headers
if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALSG')
    % Determining DataPoints
    f.BOData = f.EOexH;
    f.EOData = f.EOF;
    NSx.MetaTags.DataPoints = (f.EOF-f.EOexH)/(ChannelCount*2);
elseif or(strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD'), strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP'))
    segmentCount = 0;
    while double(ftell(FID)) < f.EOF
        if (fread(FID, 1, 'uint8') ~= 1)
            % Fixing another bug in Central 6.01.00.00 TOC where DataPoints is
            % not written back into the Data Header
            %% BIG NEEDS TO BE FIXED
            NSx.MetaTags.DataPoints = double(f.EOF - f.BOData)/(ChannelCount*2);
            break;
        end
        segmentCount = segmentCount + 1;
        %%% MODIFY THIS LINE BELOW %%%
        if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD')
            startTimeStamp = fread(FID, 1, 'uint32');
        elseif strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP')
            startTimeStamp = fread(FID, 1, 'uint64');
        end
        if strcmpi(multinsp, 'yes')
            startTimeStamp = startTimeStamp + syncShift;
            fseek(FID, -timeStampBytes, 'cof');
            fwrite(FID, startTimeStamp, '*uint32');
        end
        NSx.MetaTags.Timestamp(segmentCount)  = startTimeStamp;
        NSx.MetaTags.DataPoints(segmentCount) = fread(FID, 1, 'uint32');
        f.BOData(segmentCount) = double(ftell(FID));
        fseek(FID, NSx.MetaTags.DataPoints(segmentCount) * ChannelCount * 2, 'cof');
        f.EOData(segmentCount) = double(ftell(FID));
        % Fixing the bug in 6.01.00.00 TOC where DataPoints is not
        % updated and is left as 0
        % NSx.MetaTags.DataPoints(segmentCount) = (f.EOData(segmentCount)-f.BOData(segmentCount))/(ChannelCount*2);
    end
end

% Determining if the file has a pause in it
if length(NSx.MetaTags.DataPoints) > 1
    NSx.RawData.PausedFile = 1;
%     if modifiedTime == 1
%         disp('This data file contains pauses.');
%         disp('openNSx cannot read files with pauses using the ''t:XX'' parameter.');
%         fclose(FID); clear variables; if nargout; varargout{1} = -1; end; return;
%     end
end


%% Added by NH - Feb 19, 2014
% Create incrementing loop to skip from dataheader to dataheader and 
% collect the dataheader data in individual cells
headerCount = 0;
if NSx.RawData.PausedFile == 1
    while double(ftell(FID)) < f.EOF
        headerCount = headerCount + 1;
        fseek(FID, f.EOexH, 'bof');
        DataHeader{headerCount} = fread(FID, 9, '*uint8');
        DataPoints(headerCount) = typecast(DataHeader{headerCount}(6:9), 'uint32');

        f.BOData(headerCount) = double(ftell(FID));
        fseek(FID, DataPoints(headerCount) * ChannelCount * 2, 'cof');
        f.EOData(headerCount) = double(ftell(FID));
    end

    % Create an array that will contain all of the dataheader data
    % collected in the cells above
    FinalDataHeader = [];

    %Fill the above mentioned pre-created array
    for i = 1:headerCount
        FinalDataHeader = cat(1,FinalDataHeader,DataHeader(i));
    end

    % Convert to correct type for interpreting in separatingPausedNSx
    FinalDataHeader = cell2mat(FinalDataHeader);

    NSx.RawData.DataHeader = FinalDataHeader;

    fseek(FID, f.EOexH, 'bof');
end

%% Copying ChannelID to MetaTags for filespec 2.2, 2.3, and 3.0 for compatibility with filespec 2.1
if or(strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD'), strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP'))
    NSx.MetaTags.ChannelID = [NSx.ElectrodesInfo.ElectrodeID]';
end

%% Determining the number of channels to read and validating the input
if ~elecReading
    if ~exist('userRequestedChanRow', 'var')
        userRequestedChannels = NSx.MetaTags.ChannelID;
    else
        if any(userRequestedChanRow > ChannelCount)
            disp(['Channel file only contains ' num2str(ChannelCount) ' channels.']);
            fclose(FID); clear variables; if nargout; varargout{1} = -1; end; return;
        else
            userRequestedChannels = NSx.MetaTags.ChannelID(userRequestedChanRow);
            NSx.MetaTags.ChannelCount = length(userRequestedChannels);
        end
    end
else
    NSx.MetaTags.ChannelCount = length(userRequestedChannels);
end

for idx = 1:length(userRequestedChannels)
    if ~any(ismember(NSx.MetaTags.ChannelID, userRequestedChannels(idx)))
        disp(['Electrode ' num2str(Mapfile.Channel2Electrode(userRequestedChannels(idx))) ' does not exist in this file.']);
        fclose(FID); 
        clear variables; 
        if nargout; varargout{1} = -1; end
        return;
    end
    userRequestedChanRow(idx) = find(NSx.MetaTags.ChannelID == userRequestedChannels(idx),1);
end

%% Removing extra ElectrodesInfo for channels not read
if or(strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD'), strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP'))
    for headerIDX = length(NSx.ElectrodesInfo):-1:1
        if ~ismember(headerIDX, userRequestedChanRow)
            NSx.ElectrodesInfo(headerIDX) = [];
        end
    end
end

%% Adjusts StartPacket and EndPacket based on what time setting (sec, min,
%  hour, or packets) the user has indicated in the input argument.
if ~exist('EndPacket', 'var')
    EndPacket = sum(NSx.MetaTags.DataPoints);
end
switch TimeScale
    case 'sec'
        StartPacket = double(StartPacket) * NSx.MetaTags.SamplingFreq + 1;
        EndPacket = EndPacket * NSx.MetaTags.SamplingFreq;
    case 'min'
        StartPacket = StartPacket * NSx.MetaTags.SamplingFreq * 60 + 1;
        EndPacket = EndPacket * NSx.MetaTags.SamplingFreq * 60;
    case 'hour'
        StartPacket = StartPacket * NSx.MetaTags.SamplingFreq * 3600 + 1;
        EndPacket = EndPacket * NSx.MetaTags.SamplingFreq * 3600;
end

%% Validate StartPacket and EndPacket to make sure they do not exceed the
%  length of packets in the file. If EndPacket is over then the last packet
%  will be set for EndPacket. If StartPacket is over then will exist with an
%  error message.
if StartPacket >= EndPacket
    disp('The starting packet is greater than the end packet.');
    disp('The file was not read.');
    fclose(FID);
    if nargout; varargout{1} = -1; end
    return;
end
if StartPacket <= 0
    disp('The starting packet must be greater or equal to 1.');
    disp('The starting packet was changed to 1.');
    StartPacket = 1;
end
if EndPacket > sum(NSx.MetaTags.DataPoints)
    if StartPacket >= NSx.MetaTags.DataPoints
        disp('The starting packet is greater than the total data duration.');
        disp('The file was not read.');
        fclose(FID);
        if nargout; varargout{1} = -1; end
        return;
    end
    disp('The time interval specified is longer than the data duration.');
    disp('Last data point will be used instead.');
    disp('Press enter to continue...');
    pause;
    EndPacket = sum(NSx.MetaTags.DataPoints) - 1;
end

% Adjusting the endPacket for the skipFactor to reduce the length of
% the data read.
% DEBUG: This is not needed since the same length of data is to be
% read.
EndPacket = EndPacket / skipFactor; 

% Finding which data segment the StartPacket is falling in-between
segmentCounters = [];
startTimeStampShift = 0;
if NSx.RawData.PausedFile
    dataPointOfInterest = StartPacket;
    for idx = 1:length(NSx.MetaTags.DataPoints)
        if dataPointOfInterest <= sum(NSx.MetaTags.DataPoints(1:idx)) %sum(NSx.MetaTags.DataPoints(1:idx)) < dataPointOfInterest && ...
            if isempty(segmentCounters)
                if idx == 1
                    segmentStartPacket(idx) = dataPointOfInterest;
                else
                    segmentStartPacket(idx) = dataPointOfInterest - (sum(NSx.MetaTags.DataPoints(1:idx-1)));
                end
                startTimeStampShift = segmentStartPacket(idx) - 1;
                if EndPacket <= sum(NSx.MetaTags.DataPoints(1:idx))
                    segmentDataPoints(idx) = EndPacket - sum(NSx.MetaTags.DataPoints(1:idx-1)) - segmentStartPacket(idx) + 1;
                    segmentCounters = [idx, idx];
                    break;
                end
                segmentDataPoints(idx) = sum(NSx.MetaTags.DataPoints(1:idx)) - dataPointOfInterest + 1;
                dataPointOfInterest = EndPacket;
            else
                segmentStartPacket(idx) = 1;
                if idx == 1
                	segmentDataPoints(idx) = dataPointOfInterest;
                else
                    segmentDataPoints(idx) = dataPointOfInterest - sum(NSx.MetaTags.DataPoints(1:idx-1));
                    segmentCounters(length(segmentCounters)+1) = idx;
                end
                break;
            end
            segmentCounters(length(segmentCounters)+1) = idx;
        else
            if isempty(segmentCounters)
                segmentStartPacket(idx) = NSx.MetaTags.DataPoints(idx);
                segmentDataPoints(idx) = 0;
            elseif length(segmentCounters) == 1
                segmentStartPacket(idx) = 1;
                segmentDataPoints(idx) = NSx.MetaTags.DataPoints(idx);
            else
                segmentStartPacket(idx) = 1;
                segmentDataPoints(idx) = 0;
            end
        end
    end
end

DataLength = EndPacket - StartPacket + 1;
    
% from now StartPacket and EndPacket are in terms of Samples and are zero-based
clear TimeScale

%% Reading the data if flag 'read' is used
if strcmp(ReadData, 'read')

    % Determine what channels to read
    numChansToRead = double(length(min(userRequestedChanRow):max(userRequestedChanRow)));
    if NSx.RawData.PausedFile
        for dataIDX = segmentCounters(1):segmentCounters(2) %1:length(NSx.MetaTags.DataPoints)
            fseek(FID, f.BOData(dataIDX), 'bof');
            % Skip the file to the beginning of the time requsted, if not 0
            fseek(FID, (segmentStartPacket(dataIDX) - 1) * 2 * ChannelCount, 'cof');
            % Skip the file to the first channel to read
            fseek(FID, (find(NSx.MetaTags.ChannelID == min(userRequestedChannels))-1) * 2, 'cof');        
            % Read data
            NSx.Data{size(NSx.Data,2)+1} = fread(FID, [numChansToRead segmentDataPoints(dataIDX)], [num2str(numChansToRead) precisionType], double((ChannelCount-numChansToRead)*2 + ChannelCount*(skipFactor-1)*2));
        end
        NSx.MetaTags.DataPoints = segmentDataPoints(segmentCounters(1):segmentCounters(2));
        NSx.MetaTags.Timestamp = NSx.MetaTags.Timestamp(segmentCounters(1):segmentCounters(2));
        NSx.MetaTags.Timestamp(1) = NSx.MetaTags.Timestamp(1) + startTimeStampShift;
    else
        fseek(FID, f.BOData(1), 'bof');
        % Skip the file to the beginning of the time requsted, if not 0
        fseek(FID, (StartPacket - 1) * 2 * ChannelCount, 'cof');
        % Skip the file to the first channel to read
        fseek(FID, (find(NSx.MetaTags.ChannelID == min(userRequestedChannels))-1) * 2, 'cof');        
        % Read data
        NSx.Data = fread(FID, [numChansToRead DataLength], [num2str(numChansToRead) precisionType], double((ChannelCount-numChansToRead)*2 + ChannelCount*(skipFactor-1)*2));
    end
end

%% Fixing a bug in 6.03 TOC where an extra 0-length packet is introduced
if NSx.RawData.PausedFile && strcmp(ReadData, 'read')
    if isempty(NSx.Data{1})
        NSx.Data = cell2mat(NSx.Data(2));
    end
end

% Fixing a bug in 6.03 where data packets with 0 lengh may be added
if any(NSx.MetaTags.DataPoints == 0) && strcmp(ReadData, 'read')
    segmentsThatAreZero = find(NSx.MetaTags.DataPoints == 0);
    NSx.MetaTags.DataPoints(segmentsThatAreZero) = [];
    NSx.MetaTags.Timestamp(segmentsThatAreZero) = [];
    NSx.Data(segmentsThatAreZero) = [];
end

%% Removing extra channels that were read, but weren't supposed to be read
channelThatWereRead = min(userRequestedChanRow):max(userRequestedChanRow);
if ~isempty(setdiff(channelThatWereRead,userRequestedChanRow))
	deleteChannels = setdiff(channelThatWereRead, userRequestedChanRow) - min(userRequestedChanRow) + 1;
    if NSx.RawData.PausedFile
        for segIDX = 1:size(NSx.Data,2)
            NSx.Data{segIDX}(deleteChannels,:) = [];
        end
    else
        NSx.Data(deleteChannels,:) = [];
    end
end

%% Remove the cell if there is only one recorded segment present
if all(size(NSx.Data) == [1,1])
    NSx.Data = NSx.Data{1};
end

%% Adjusting the ChannelID variable to match the read electrodes
channelIDToDelete = setdiff(1:ChannelCount, userRequestedChanRow);
NSx.MetaTags.ChannelID(channelIDToDelete) = [];

%% Adjusting the file for a non-0 timestamp start
if ~NSx.RawData.PausedFile & StartPacket == 1
    if length(NSx.MetaTags.Timestamp) > 1
        cellIDX = 1; % only do this for the first cell segment and not modify the subsequent segments
        if strcmpi(ReadData, 'read')
            NSx.Data{cellIDX} = [zeros(NSx.MetaTags.ChannelCount, floor(NSx.MetaTags.Timestamp(cellIDX) / skipFactor)) NSx.Data{cellIDX}];
        end
        NSx.MetaTags.DataPoints(cellIDX) = NSx.MetaTags.DataPoints(cellIDX) + NSx.MetaTags.Timestamp(cellIDX);
        NSx.MetaTags.DataDurationSec(cellIDX) = NSx.MetaTags.DataPoints(cellIDX) / NSx.MetaTags.SamplingFreq;
        NSx.MetaTags.Timestamp(cellIDX) = 0;
    elseif strcmpi(ReadData, 'read')
        NSx.Data = [zeros(NSx.MetaTags.ChannelCount, floor(NSx.MetaTags.Timestamp / skipFactor)) NSx.Data];
        NSx.MetaTags.DataPoints = size(NSx.Data,2);
        NSx.MetaTags.DataDurationSec = NSx.MetaTags.DataPoints / NSx.MetaTags.SamplingFreq;
        NSx.MetaTags.Timestamp = 0;
    end

    if strcmpi(multinsp, 'yes')
        NSx.Data = [zeros(NSx.MetaTags.ChannelCount, syncShift) NSx.Data];
        NSx.MetaTags.DataPoints = size(NSx.Data,2);
        NSx.MetaTags.DataDurationSec = NSx.MetaTags.DataPoints / NSx.MetaTags.SamplingFreq;
    end
end

%% Adjusting for the data's unit.
if strcmpi(waveformUnits, 'uV')
    if iscell(NSx.Data) % Contribution by Michele Cox @ Vanderbilt
    	NSx.Data = cellfun(@(x) bsxfun(@rdivide, double(x), 1./(double([NSx.ElectrodesInfo.MaxAnalogValue])./double([NSx.ElectrodesInfo.MaxDigiValue]))'),NSx.Data ,'UniformOutput',false);
    else
        NSx.Data = bsxfun(@rdivide, double(NSx.Data), 1./(double([NSx.ElectrodesInfo.MaxAnalogValue])./double([NSx.ElectrodesInfo.MaxDigiValue]))');
    end % End of contribution
end

%% Converting the data points in sample to seconds
NSx.MetaTags.DataPointsSec = double(NSx.MetaTags.DataPoints)/NSx.MetaTags.SamplingFreq;

%% Calculating the DataPoints in seconds and adding it to MetaData
NSx.MetaTags.DataDurationSec = double(NSx.MetaTags.DataPoints)/NSx.MetaTags.SamplingFreq;

%% Displaying a report of basic file information and the Basic Header.
if strcmp(Report, 'report')
    disp( '*** FILE INFO **************************');
    disp(['File Path          = '  NSx.MetaTags.FilePath]);
    disp(['File Name          = '  NSx.MetaTags.Filename]);
    disp(['File Extension     = '  NSx.MetaTags.FileExt]);
	disp(['File Version       = '  NSx.MetaTags.FileSpec]);
    disp(['Duration (seconds) = '  num2str(NSx.MetaTags.DataDurationSec)]);
    disp(['Total Data Points  = '  num2str(NSx.MetaTags.DataPoints)]);
    disp(' ');
    disp( '*** BASIC HEADER ***********************');
    disp(['File Type ID       = '  NSx.MetaTags.FileTypeID]);
    disp(['Sample Frequency   = '  num2str(double(NSx.MetaTags.SamplingFreq))]);
    disp(['Electrodes Read    = '  num2str(double(NSx.MetaTags.ChannelCount))]);
    disp(['Data Point Read    = '  num2str(size(NSx.Data,2))]);
end

%% If user does not specify an output argument it will automatically create
%  a structure.
outputName = ['NS' fext(4)];
if (nargout == 0)
    assignin('caller', outputName, NSx);
else
    varargout{1} = NSx;
end

if strcmp(Report, 'report')
	disp(['The load time for ' outputName ' file was ' num2str(toc, '%0.1f') ' seconds.']);
end
fclose(FID);

end