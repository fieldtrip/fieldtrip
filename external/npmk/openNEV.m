function varargout = openNEV(varargin)

% openNEV
%
% Opens an .nev file for reading, returns all file information in a NEV
% structure. Works with File Spec 2.1 & 2.2 & 2.3 & 3.0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Use OUTPUT = openNEV(fname, 'noread', 'report', 'noparse', 'nowarning', 
%                             'nosave', 'nomat', 'uV', 'overwrite', 'direct').
% 
% NOTE: All input arguments are optional. Input arguments may be in any order.
%
%   fname:        Name of the file to be opened. If the fname is omitted
%                 the user will be prompted to select a file using an open 
%                 file user interface. 
%                 DEFAULT: Will open Open File UI.
%
%   'noread':     Will not read the spike waveforms if user passes this argument.
%                 DEFAULT: will read spike waveform.
%
%   'report':     Will show a summary report if user passes this argument.
%                 DEFAULT: will not show report.
%
%   'parse':    The code will not parse the experimental parameters in digital I/O.
%                 See below for guidelines on how to format your parameters.
%                 DEFAULT: will not parse the parameters.
%
%   'nowarning':  The code will not give a warning if there is an error in
%                 parsing.
%                 DEFAULT: will give warning message.
%
%   'nosave':     The code will not save a copy of the NEV structure as a
%                 MAT file. By default the code will save a copy in the same
%                 folder as the NEV file for easy future access.
%                 DEFAULT: will save the MAT file.
%
%   'nomat':      Will not look for a MAT file. This option will force
%                 openNEV to open a NEV file instead of any available MAT
%                 files.
%                 DEFAULT: will load the MAT file if available.
%
%   'uV':         Will read the spike waveforms in unit of uV instead of
%                 raw values. Note that this conversion may lead to loss of
%                 information (e.g. 15/4 = 4) since the waveforms type will
%                 stay in int16. It's recommended to read raw spike
%                 waveforms and then perform the conversion at a later
%                 time.
%                 DEFAULT: will read waveform information in raw.
%
%   '8bits':      Indicates that 8 bits on the digital IO port was used
%                 instead of 16 bits.
%                 DEFAULT: will assumes that 16 bits of digital IO were used.
%
%   't:':         Indicats the time window of the NEV file to be read. For
%                 example, if t: is set to 2 (i.e. 't:0.6')
%                 then only the first 2 seconds of the file is to be read. If set
%                 to 2-50 (i.e. 't:2:50) then the time between 2 seconds
%                 and 50 seconds will be read.
%                 DEFAULT: the entire file will be read if 't:xx:xx' is not
%                 passed to the function.
%
%   'overwrite':  If MATLAB loads a NEV file using 'nomat' and a MAT file
%                 already exists, by default it will prompt the user to
%                 allow for overwriting the old MAT. Passing the
%                 'overwrite' flag will automatically overwrite the newly
%                 opened NEV file ont the old MAT file.
%                 DEFAULT: will ask the user whether to overwrite the old
%                 MAT.
%
%   'direct':     Use this if you are using a CerePlex Direct system
%                 without the typical strobe mode. This will treat the 16th
%                 bit of the digital input as a strobe signal and report
%                 the remaining 15 bits as the digital input value. 
%
%   OUTPUT:       Contains the NEV structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   USAGE EXAMPLE: 
%   
%   openNEV('report','read');
%
%   In the example above, the file dialogue will prompt for a file. A
%   report of the file contents will be shown. The digital data will not be
%   parsed. The data needs to be in the proper format (refer below). The 
%   spike waveforms are in raw units and not in uV.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DIGITAL PARAMETERS/MARKERS FORMAT:
%
%   In order for this function to parse your experimental parameters they 
%   need to be in the following format:
%
%   *ParamLabel:Parameter1=value1;Parameter2=value2;Parameter3=value3;#
%
%   TWO EXAMPLES:
%   *ExpParameter:Intensity=1.02;Duration=400;Trials=1;PageSegment=14;#
%
%   *Stimulation:StimCount=5;Duration=10;#
%
%   In the first example, the parameter is of type "ExpParameter". The 
%   parameters are, "Intensity, Duration, Trials, and PageSement." The 
%   values of those parameters are, "1.02, 400, 1, and 14," respectively.
%   The second example is of type "Stimulation". The name of the parameters
%   are "StimCount" and "Duration" and the values are "5" and "10" 
%   respectively.
%   -----------------------------------------------------------------------
%   It can also read single value markers that follow the following format.
%
%   *MarkerName=Value;#
%
%   EXAMPLES:  *WaitSeconds=10;# OR  *JuiceStatus=ON;#
%
%   The above line is a "Marker". The marker value is 10 in the first 
%   and it's ON in the second example.
%   -----------------------------------------------------------------------
%   Moreover, the marker could be a single value:
%
%   *MarkerValue#
%
%   EXAMPLES: *JuiceOff#  OR  *HandsOnSwitches#
%   -----------------------------------------------------------------------
%   The label, parameter name, and values are flexible and can be anything.
%   The only required formatting is that the user needs to have a label
%   followed by a colon ':', followed by a field name 'MarkerVal', followed
%   by an equal sign '=', followed by the parameter value '10', and end
%   with a semi-colon ';'.
%
%   NOTE:
%   Every line requires a pound-sign '#' at the very end. 
%   Every line requires a star sign '*' at the very beginning. If you
%   use LabVIEW SendtoCerebus.vi by Kian Torab then there is no need for 
%   a '*' in the beginning.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 4.4.0.0:
%   - Major performance boost in reading NEV files when tracking data is 
%   stored in the file.
%
% 4.4.0.2:
%   - Updated documentation.
%
% 4.4.0.3: 5 January 2014
%   - Fixed the way DayOfWeek is read in MetaTags.
%   - Fixed 'noread' argument, so when passed, openNEV will not read the
%   spike waveforms.
%
% 4.4.1.0: 25 January 2014
%   - Fixed a bug that resulted from passing 'read' to openNEV.
%
% 4.4.2.0: 28 February 2014
%   - Fixed bug related to loading data with t:XX:XX argument.
%
% 4.4.3.0: 12 June 2014
%   - Fixed a typo in the help.
%
% 4.4.3.1: 13 June 2014
%   - Updated the version numbers in the help and in the function itself.
%
% 5.0.0.0: 02 December 2014
%   - Fixed a bug where Application name wasn't being read properly.
%   - Warnings now don't show up in more places when "nowarning" is used.
%   - Added field FileExt to MetaTags.
%   - Added 512 synchronized reading capability.
%   - Fixed the date in NSx.MetaTags.DateTime.
%
% 5.1.0.0: 28 March 2015
%   - Added the ability to read from networked drives in Windows.
%   - Fixed the DateTime variable in MetaTags.
%   - Fixed the date in NSx.MetaTags.DateTime (again).
%   - Fixed a bug related to >512-ch data loading.
%
% 5.1.1.0: 1 April 2015
%   - Fixed a bug with NeuroMotive when spike window is changed from the
%     original length.
%
% 5.1.2.0: June 30 2015
%   - Fixed a bug regarding the number of packages when 'no read' is used.
%
% 5.1.3.0: July 10 2015
%   - Fixed a bug with NeuroMotive data reading when both objects and
%     markers were being recorded.
%
% 5.2.0.0: June 11 2016
%   - Added support for CerePlex Direct strobe mode on digital input.
%   - Fixed a bug with reading NeuroMotive data that resulted in a crash.
%
% 5.3.0.0: June 13 2016
%   - Fixed a bug with reading NeuroMotive data that resulted in a crash.
%   - Improved and more detailed parsing of NeuroMotive events.
%   - Added parsing of comment start time and comment committ time (time
%     that a comment is entered.
%
% 5.3.1.0: September 1, 2017
%   - Fixed a bug with file path and whent this was passed to the function.
%
% 5.4.0.0: September 13, 2017
%   - Checks to see if there's a newer version of NPMK is available.
%   - Properly reads the comment colors.
%
% 5.4.0.1: January 10, 2018
%   - Fixed a NeuroMotive bug when AllMarkers was being recorded.
%
% 5.4.1.0: April 25, 2018
%   - Now all comments open in order.
%   - Fixed a bug with path of file if both NEV and MAT were moved to a new
%     location.
%
% 6.0.0.0: January 27, 2020
%   - Added support for 64-bit timestamps in NEV and NSx.
%   - Removed dependency on MATLAB R2016b by removing function 'contains'.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check for the latest version fo NPMK
NPMKverChecker

%% Defining structures
NEV = struct('MetaTags',[], 'ElectrodesInfo', [], 'Data', []);
NEV.MetaTags = struct('Subject', [], 'Experimenter', [], 'DateTime', [],...
    'SampleRes',[],'Comment',[],'FileTypeID',[],'Flags',[], 'openNEVver', [], ...
    'DateTimeRaw', [], 'FileSpec', [], 'PacketBytes', [], 'HeaderOffset', [], ...
    'DataDuration', [], 'DataDurationSec', [], 'PacketCount', [], ...
    'TimeRes', [], 'Application', [], 'Filename', [], 'FilePath', []);
NEV.Data = struct('SerialDigitalIO', [], 'Spikes', [], 'Comments', [], 'VideoSync', [], ...
    'Tracking', [], 'TrackingEvents', [], 'PatientTrigger', [], 'Reconfig', []);
NEV.Data.Spikes = struct('TimeStamp', [],'Electrode', [],...
    'Unit', [],'Waveform', [], 'WaveformUnit', []);
NEV.Data.SerialDigitalIO = struct('InputType', [], 'TimeStamp', [],...
    'TimeStampSec', [], 'Type', [], 'Value', [], 'InsertionReason', [], 'UnparsedData', []);
NEV.Data.VideoSync = struct('TimeStamp', [], 'FileNumber', [], 'FrameNumber', [], 'ElapsedTime', [], 'SourceID', []);
NEV.Data.Comments = struct('TimeStampStarted', [], 'TimeStampStartedSec', [], 'TimeStamp', [], 'TimeStampSec', [], 'CharSet', [], 'Text', []);
NEV.Data.Tracking = [];
NEV.Data.TrackingEvents = struct('TimeStamp', [], 'TimeStampSec', [], 'ROIName', [], 'ROINum', [], 'Event', [], 'Frame', []);
NEV.Data.PatientTrigger = struct('TimeStamp', [], 'TriggerType', []);
NEV.Data.Reconfig = struct('TimeStamp', [], 'ChangeType', [], 'CompName', [], 'ConfigChanged', []);
Flags = struct;

%% Check for multiple versions of openNEV in path
if size(which('openNEV', '-ALL'),1) > 1
    disp('WARNING: There are multiple openNEV functions in the path. Use which openNEV -ALL for more information.');
end

%% Validating input arguments
for i=1:length(varargin)
    switch lower(varargin{i})
        case 'report'
            Flags.Report = varargin{i};
        case 'noread'
            Flags.ReadData = varargin{i};
        case 'nomultinsp'
            Flags.MultiNSP = varargin{i};
        case 'read'
            Flags.ReadData = varargin{i};
        case 'nosave'
            Flags.SaveFile = varargin{i};
        case 'nomat'
            Flags.NoMAT = varargin{i};
        case 'direct'
            Flags.Direct = varargin{i};
        case 'nowarning'
            Flags.WarningStat = varargin{i};
        case 'parse'
            Flags.ParseData = 'parse';
        case 'uv'
            Flags.waveformUnits = 'uV';
        case '8bits'
            Flags.digIOBits = '8bits';
        case 'overwrite'
            Flags.Overwrite = 'overwrite';
        case 'nooverwrite'
            Flags.Overwrite = 'nooverwrite';
        otherwise
            temp = varargin{i};
            if length(temp)>3 && ...
                    (strcmpi(temp(3),'\') || ...
                     strcmpi(temp(1),'/') || ...
                     strcmpi(temp(2),'/') || ...
                     strcmpi(temp(1:2), '\\') || ...
                     strcmpi(temp(end-3), '.'))
                fileFullPath = varargin{i};
                if exist(fileFullPath, 'file') ~= 2
                    disp('The file does not exist.');
                    varargout{1} = [];
                    return;
                end
            elseif length(temp)>3 && strcmpi(temp(1:2),'t:') && ~strcmpi(temp(3), '\') && ~strcmpi(temp(3), '/')
                temp(1:2) = [];
                temp = str2num(temp);
                if length(temp) == 1
                    fprintf('Only one timepoint (%0.0f) was passed to the function.\n', temp);
                    fprintf('The initial timepoint is set to 0, so data between 0 and %0.0f will be read.\n', temp);
                    temp(2) = temp;
                    temp(1) = 0;
                end
                readTime = [temp(1), temp(end)];
                Flags.SaveFile = 'nosave';
                Flags.NoMAT = 'nomat';
            elseif (strncmp(temp, 'c:', 2) && temp(3) ~= '\' && temp(3) ~= '/')
                Flags.selChannels = str2num(temp(3:end)); %#ok<ST2NM>
            else
                if ~isnumeric(varargin{i})
                    disp(['Invalid argument ''' varargin{i} ''' .']);
                else
                    disp(['Invalid argument ''' num2str(varargin{i}) ''' .']);
                end
                clear variables;
                if nargout
                    varargout{1} = [];
                end
                return;
            end
            clear temp;
    end
end; clear i;

%% Defining and validating variables
if ~exist('fileFullPath', 'var')
    if exist('getFile.m', 'file') == 2
        [fileName pathName] = getFile('*.nev;*.nevm', 'Choose a NEV file...');
    else
        [fileName pathName] = uigetfile;
    end
    fileFullPath = [pathName fileName];
    if fileFullPath==0; 
        clear variables; 
        if nargout
            varargout{1} = [];
        end
        disp('No file was selected.');
        return
    end
end

[~, ~, fileExt] = fileparts(fileFullPath);

%% Loading .x files for multiNSP configuration
if strcmpi(fileExt(2:4), 'nev') && length(fileExt) == 5
    fileFullPath(1) = fileFullPath(end);
    fileFullPath(end) = [];
end


if ~isfield(Flags, 'Report');        Flags.Report = 'noreport'; end
if ~isfield(Flags, 'WarningStat');   Flags.WarningStat = 'warning'; end;
if ~isfield(Flags, 'ReadData');      Flags.ReadData = 'read'; end
if ~isfield(Flags, 'ParseData');     Flags.ParseData = 'noparse'; end
if ~isfield(Flags, 'SaveFile');      Flags.SaveFile = 'save'; end;
if ~isfield(Flags, 'NoMAT');         Flags.NoMAT = 'yesmat'; end;
if ~isfield(Flags, 'waveformUnits'); Flags.waveformUnits = 'raw'; end;
if ~isfield(Flags, 'digIOBits');     Flags.digIOBits = '16bits'; end;
if ~isfield(Flags, 'Overwrite');     Flags.Overwrite = 'nooverwrite'; end;
if ~isfield(Flags, 'MultiNSP');      Flags.MultiNSP = 'multinsp'; end;
if ~isfield(Flags, 'selChannels');   Flags.selChannels = 'all'; end;
if ~isfield(Flags, 'Direct');        Flags.Direct = 'nodirect'; end;

if strcmpi(Flags.Report, 'report')
    disp(['openNEV ' NEV.MetaTags.openNEVver]);
end

syncShift = 0;

% Check to see if 512 setup and calculate offset
if strcmpi(Flags.MultiNSP, 'multinsp')
    fiveTwelveFlag = regexp(fileFullPath, '-i[0123]-');
    if ~isempty(fiveTwelveFlag)
        syncShift = multiNSPSync(fileFullPath);
    else
        Flags.MultiNSP = 'no';
    end
end

%%  Validating existance of parseCommand
if strcmpi(Flags.ParseData, 'parse') 
    if exist('parseCommand.m', 'file') ~= 2
        disp('This version of openNEV requires function parseCommand.m to be placed in path.');
        clear variables;
        if nargout
            varargout{1} = [];
        end
        return;
    end
end

tic;
matPath = [fileFullPath(1:end-4) '.mat'];

%% Check for a MAT file and load that instead of NEV
if exist(matPath, 'file') == 2 && strcmpi(Flags.NoMAT, 'yesmat') && strcmpi(Flags.WarningStat, 'warning')
    disp('MAT file corresponding to selected NEV file already exists. Loading MAT instead...');
    load(matPath);
    NEV.MetaTags.FilePath = fileFullPath;    
    if isempty(NEV.Data.Spikes.Waveform) && strcmpi(Flags.ReadData, 'read') && strcmpi(Flags.WarningStat, 'warning')
        disp('The MAT file does not waveforms. Loading NEV instead...');
    else
        NEV = killUnwantedChannels(NEV, Flags.selChannels);
        if ~nargout
            assignin('base', 'NEV', NEV);
            clear variables;
        else
            varargout{1} = NEV;
        end
        return;
    end
end

%% Reading BasicHeader information from file
FID                       = fopen(fileFullPath, 'r', 'ieee-le');
BasicHeader               = fread(FID, 336, '*uint8');
NEV.MetaTags.FileTypeID   = char(BasicHeader(1:8)');
NEV.MetaTags.FileSpec     = [num2str(double(BasicHeader(9))) '.' num2str(double(BasicHeader(10)))];
NEV.MetaTags.Flags        = dec2bin(double(typecast(BasicHeader(11:12), 'uint16')),16);
Trackers.fExtendedHeader  = double(typecast(BasicHeader(13:16), 'uint32'));
NEV.MetaTags.HeaderOffset = Trackers.fExtendedHeader;
Trackers.countPacketBytes = double(typecast(BasicHeader(17:20), 'uint32'));
NEV.MetaTags.PacketBytes  = Trackers.countPacketBytes;
NEV.MetaTags.TimeRes      = double(typecast(BasicHeader(21:24), 'uint32'));
NEV.MetaTags.SampleRes    = typecast(BasicHeader(25:28), 'uint32');
t                         = double(typecast(BasicHeader(29:44), 'uint16'));
tempApp                   = BasicHeader(45:76)';
tempApp(find(tempApp == 0):end) = [];
NEV.MetaTags.Application  = char(tempApp); clear tempApp;
NEV.MetaTags.Comment      = char(BasicHeader(77:332)');
[NEV.MetaTags.FilePath, NEV.MetaTags.Filename, NEV.MetaTags.FileExt] = fileparts(fileFullPath);
Trackers.countExtHeader   = typecast(BasicHeader(333:336), 'uint32');
clear BasicHeader;

if or(strcmpi(NEV.MetaTags.FileTypeID, 'NEURALEV'), strcmpi(NEV.MetaTags.FileTypeID, 'BREVENTS'))
    if exist([fileFullPath(1:end-8) '.sif'], 'file') == 2
        METATAGS = textread([fileFullPath(1:end-8) '.sif'], '%s');
        NEV.MetaTags.Subject      = METATAGS{3}(5:end-5);
        NEV.MetaTags.Experimenter = [METATAGS{5}(8:end-8) ' ' METATAGS{6}(7:end-7)];
    end
end
if ~any(strcmpi(NEV.MetaTags.FileSpec, {'2.1', '2.2', '2.3', '3.0'}))
    disp('Unknown filespec. Cannot open file.');
    fclose FID;
    clear variables;
    if nargout
        varargout{1} = [];
    end
    return;
end
clear fileFullPath;

%% Parsing and validating FileSpec and DateTime variables
NEV.MetaTags.DateTimeRaw = t.';
NEV.MetaTags.DateTime = datestr(datenum(t(1), t(2), t(4), t(5), t(6), t(7)));
clear t;

%% Removing extra garbage characters from the Comment field.
NEV.MetaTags.Comment(find(NEV.MetaTags.Comment==0,1):end) = 0;

%% Recording after BasicHeader file position
Trackers.fBasicHeader = ftell(FID); %#ok<NASGU>

% Calculating the length of the data
currentLocation = ftell(FID);
fseek(FID, -Trackers.countPacketBytes, 'eof');
NEV.MetaTags.DataDuration = fread(FID, 1, 'uint32=>double');
NEV.MetaTags.DataDurationSec = double(NEV.MetaTags.DataDuration) / double(NEV.MetaTags.SampleRes);
fseek(FID, currentLocation, 'bof');

%% Reading ExtendedHeader information
for ii=1:Trackers.countExtHeader
    ExtendedHeader = fread(FID, 32, '*uint8');
    PacketID = char(ExtendedHeader(1:8)');
    switch PacketID
        case 'ARRAYNME'
            NEV.ArrayInfo.ElectrodeName    = char(ExtendedHeader(9:end));
        case 'ECOMMENT'
            NEV.ArrayInfo.ArrayComment     = char(ExtendedHeader(9:end));
        case 'CCOMMENT'
            NEV.ArrayInfo.ArrayCommentCont = char(ExtendedHeader(9:end));
        case 'MAPFILE'
            NEV.ArrayInfo.MapFile          = char(ExtendedHeader(9:end));
        case 'NEUEVWAV'
            ElectrodeID                       = typecast(ExtendedHeader(9:10), 'uint16');
            NEV.ElectrodesInfo(ElectrodeID).ElectrodeID     = ElectrodeID;
            NEV.ElectrodesInfo(ElectrodeID).ConnectorBank   = char(ExtendedHeader(11)+64);
            NEV.ElectrodesInfo(ElectrodeID).ConnectorPin    = ExtendedHeader(12);
            df   = typecast(ExtendedHeader(13:14),'int16');
            % This is a workaround for the DigitalFactor overflow in NEV 
            % files. Remove once Central is updated
            if df == 21516
                NEV.ElectrodesInfo(ElectrodeID).DigitalFactor = 152592.547;
            else
                NEV.ElectrodesInfo(ElectrodeID).DigitalFactor = df;
            end
            % End of workaround
            NEV.ElectrodesInfo(ElectrodeID).EnergyThreshold = typecast(ExtendedHeader(15:16),'uint16');
            NEV.ElectrodesInfo(ElectrodeID).HighThreshold   = typecast(ExtendedHeader(17:18),'int16');
            NEV.ElectrodesInfo(ElectrodeID).LowThreshold    = typecast(ExtendedHeader(19:20),'int16');
            NEV.ElectrodesInfo(ElectrodeID).Units           = ExtendedHeader(21);
            NEV.ElectrodesInfo(ElectrodeID).WaveformBytes   = ExtendedHeader(22);
            clear ElectrodeID;
        case 'NEUEVLBL'
            ElectrodeID                       = typecast(ExtendedHeader(9:10), 'uint16');
            NEV.ElectrodesInfo(ElectrodeID).ElectrodeLabel = char(ExtendedHeader(11:26));
            clear ElectrodeID;
        case 'NEUEVFLT'
            ElectrodeID                       = typecast(ExtendedHeader(9:10), 'uint16');
            NEV.ElectrodesInfo(ElectrodeID).HighFreqCorner = typecast(ExtendedHeader(11:14),'uint32');
            NEV.ElectrodesInfo(ElectrodeID).HighFreqOrder  = typecast(ExtendedHeader(15:18),'uint32');
            NEV.ElectrodesInfo(ElectrodeID).HighFilterType = typecast(ExtendedHeader(19:20),'uint16');
            NEV.ElectrodesInfo(ElectrodeID).LowFreqCorner  = typecast(ExtendedHeader(21:24),'uint32');
            NEV.ElectrodesInfo(ElectrodeID).LowFreqOrder   = typecast(ExtendedHeader(25:28),'uint32');
            NEV.ElectrodesInfo(ElectrodeID).LowFilterType  = typecast(ExtendedHeader(29:30),'uint16');
            clear ElectrodeID;
        case 'DIGLABEL'
            Mode                    = ExtendedHeader(25);
            NEV.IOLabels{Mode+1}    = char(ExtendedHeader(9:24).');
            clear Mode;
        case 'NSASEXEV' %% Not implemented in the Cerebus firmware. 
                        %% Needs to be updated once implemented into the 
                        %% firmware by Blackrock Microsystems.
            NEV.NSAS.Freq          = typecast(ExtendedHeader(9:10),'uint16');
            NEV.NSAS.DigInputConf  = char(ExtendedHeader(11));
            NEV.NSAS.AnalCh1Conf   = char(ExtendedHeader(12));
            NEV.NSAS.AnalCh1Detect = typecast(ExtendedHeader(13:14),'uint16');
            NEV.NSAS.AnalCh2Conf   = char(ExtendedHeader(15));
            NEV.NSAS.AnalCh2Detect = typecast(ExtendedHeader(16:17),'uint16');
            NEV.NSAS.AnalCh3Conf   = char(ExtendedHeader(18));
            NEV.NSAS.AnalCh3Detect = typecast(ExtendedHeader(19:20),'uint16');
            NEV.NSAS.AnalCh4Conf   = char(ExtendedHeader(21));
            NEV.NSAS.AnalCh4Detect = typecast(ExtendedHeader(22:23),'uint16');
            NEV.NSAS.AnalCh5Conf   = char(ExtendedHeader(24));
            NEV.NSAS.AnalCh5Detect = typecast(ExtendedHeader(25:26),'uint16');
        case 'VIDEOSYN'
            cnt = 1;
            if (isfield(NEV, 'VideoSyncInfo'))
                cnt = size(NEV.VideoSyncInfo, 2) + 1;
            end
            NEV.VideoSyncInfo(cnt).SourceID     = typecast(ExtendedHeader(9:10),'uint16');
            NEV.VideoSyncInfo(cnt).SourceName   = char(ExtendedHeader(11:26))';
            NEV.VideoSyncInfo(cnt).FrameRateFPS = typecast(ExtendedHeader(27:30),'single')';
            clear cnt;
        case 'TRACKOBJ'
            cnt = 1;
            if (isfield(NEV, 'ObjTrackInfo'))
                cnt = size(NEV.ObjTrackInfo, 2) + 1;
            end
            NEV.ObjTrackInfo(cnt).TrackableType = typecast(ExtendedHeader(9:10),'uint16');
            NEV.ObjTrackInfo(cnt).TrackableID   = typecast(ExtendedHeader(11:14), 'uint32');
            NEV.ObjTrackInfo(cnt).TrackableName = char(ExtendedHeader(15:30))';
            clear cnt;
        otherwise
            disp(['PacketID ' PacketID ' is invalid.']);
            disp('Please make sure this version of openNEV is compatible with your current NSP firmware.')
            fclose(FID);
            clear variables; 
            if nargout
                varargout{1} = [];
            end
            return;
    end
end
NEV.MetaTags.ChannelID = [NEV.ElectrodesInfo.ElectrodeID];
clear ExtendedHeader PacketID ii;

%% Recording after ExtendedHeader file position and calculating Data Length
%  and number of data packets
fseek(FID, 0, 'eof');
Trackers.fData = ftell(FID);
Trackers.countDataPacket = (Trackers.fData - Trackers.fExtendedHeader)/Trackers.countPacketBytes;
NEV.MetaTags.PacketCount = Trackers.countDataPacket;

%%
Flags.UnparsedDigitalData = 0;

%% Reading packet headers and digital values
Timestamp = [];
PacketIDs = [];
tempClassOrReason = [];
tempDigiVals = [];
if NEV.MetaTags.PacketCount ~= 0
    fseek(FID, Trackers.fExtendedHeader, 'bof');
    if strcmpi(NEV.MetaTags.FileTypeID, 'NEURALEV')
        tRawData  = fread(FID, [10 Trackers.countDataPacket], '10*uint8=>uint8', Trackers.countPacketBytes - 10);
        Timestamp = tRawData(1:4,:);
        Timestamp = typecast(Timestamp(:), 'uint32').' + syncShift;
        timeStampBytes = 4;
    elseif strcmpi(NEV.MetaTags.FileTypeID, 'BREVENTS')
        tRawData  = fread(FID, [14 Trackers.countDataPacket], '14*uint8=>uint8', Trackers.countPacketBytes - 14);
        Timestamp = tRawData(1:8,:);
        Timestamp = typecast(Timestamp(:), 'uint64').' + syncShift;
        timeStampBytes = 8;
    end

    %% Calculate the number of packets that need to be read based on the time input parameters
    if ~exist('readTime', 'var')
        Trackers.readPackets = [1, length(Timestamp)];
    else
        [tmp,tempReadPackets] = find(Timestamp > readTime(1)*NEV.MetaTags.SampleRes,1,'first');
        if ~isempty(tempReadPackets)
            Trackers.readPackets(1) = tempReadPackets;
        else
            Trackers.readPackets(1) = NaN;
        end
        if isnan(Trackers.readPackets(1))
            fprintf('The file contains %0.2f seconds of data. The requested begining timestamp of %0.2f seconds is longer than the duration.\n', ...
                double(Timestamp(end))/double(NEV.MetaTags.SampleRes), ...
                readTime(1));
            clear variables;
            if nargout
                varargout{1} = [];
            end
            return;
        end
        [tmp,tempReadPackets] = find(Timestamp < readTime(2)*NEV.MetaTags.SampleRes,1,'last');
        if ~isempty(tempReadPackets)
            if readTime(2)*NEV.MetaTags.SampleRes > Timestamp(end)
                fprintf('The file contains %0.2f seconds of data. The requested end duration of %0.2f seconds will be adjusted to %0.2f seconds.\n', ...
                    double(Timestamp(end))/double(NEV.MetaTags.SampleRes), ...
                    readTime(2),...
                    double(Timestamp(end))/double(NEV.MetaTags.SampleRes));
            end
            Trackers.readPackets(2) = tempReadPackets;
        else
            Trackers.readPackets(2) = NaN;
        end
        clear tmp, tempReadPackets;
    end
    PacketIDs = tRawData(timeStampBytes+1:timeStampBytes+2,Trackers.readPackets(1):Trackers.readPackets(2));
    PacketIDs = typecast(PacketIDs(:), 'uint16').';
    tempClassOrReason = uint8(tRawData(timeStampBytes+3,Trackers.readPackets(1):Trackers.readPackets(2)));
    if strcmpi(Flags.digIOBits, '16bits')
        tempDigiVals      = tRawData(timeStampBytes+5:timeStampBytes+6,Trackers.readPackets(1):Trackers.readPackets(2));
        tempDigiVals      = typecast(tempDigiVals(:), 'uint16');
    else
        tempDigiVals      = uint16(tRawData(timeStampBytes+5,Trackers.readPackets(1):Trackers.readPackets(2)));
    end
    clear tRawData;
else
    Trackers.readPackets = zeros(1,2);
end

%% Defining PacketID constants
digserPacketID = 0;
neuralIndicesPacketIDBounds = [1, 16384];
commentPacketID = 65535;
videoSyncPacketID = 65534;
trackingPacketID = 65533;
patientTrigPacketID = 65532;
logEventPacketID = 65531;
reconfigPacketID = 65530;
recStartPacketID = 65529;


%% Parse read digital data. Please refer to help to learn about the proper
% formatting if the data.
digserIndices              = find(PacketIDs == digserPacketID);
neuralIndices              = find(neuralIndicesPacketIDBounds(2) >= PacketIDs & PacketIDs >= neuralIndicesPacketIDBounds(1));
commentIndices             = find(PacketIDs == commentPacketID);
videoSyncPacketIDIndices   = find(PacketIDs == videoSyncPacketID);
trackingPacketIDIndices    = find(PacketIDs == trackingPacketID);
patientTrigPacketIDIndices = find(PacketIDs == patientTrigPacketID);
logEventPacketIDIndices    = find(PacketIDs == logEventPacketID);
reconfigPacketIDIndices    = find(PacketIDs == reconfigPacketID);
recStartPacketIDIndices    = find(PacketIDs == recStartPacketID);
clear digserPacketID neuralIndicesPacketIDBounds commentPacketID ...
      videoSyncPacketID trackingPacketID patientTrigPacketID reconfigPacketID;
digserTimestamp            = Timestamp(digserIndices);
NEV.Data.Spikes.TimeStamp  = Timestamp(neuralIndices);
NEV.Data.Spikes.Electrode  = PacketIDs(neuralIndices);
clear PacketIDs;
NEV.Data.Spikes.Unit       = tempClassOrReason(neuralIndices); 
%clear neuralIndices;
NEV.Data.SerialDigitalIO.InsertionReason   = tempClassOrReason(digserIndices);
clear tempClassOrReason;
DigiValues                 = tempDigiVals(digserIndices);
clear tempDigiVals;

%% Reads the waveforms if 'read' is passed to the function
if strcmpi(Flags.ReadData, 'read')
    allExtraDataPacketIndices  = [commentIndices, ...
                                  videoSyncPacketIDIndices, ...
                                  trackingPacketIDIndices, ...
                                  patientTrigPacketIDIndices, ...
                                  logEventPacketIDIndices,...
                                  reconfigPacketIDIndices,...
                                  recStartPacketIDIndices];
      
    if ~isempty(allExtraDataPacketIndices) % if there is any extra packets
        fseek(FID, Trackers.fExtendedHeader, 'bof');
        fseek(FID, (Trackers.readPackets(1)-1) * Trackers.countPacketBytes, 'cof');
        tRawData  = fread(FID, [Trackers.countPacketBytes Trackers.readPackets(2)], ...
            [num2str(Trackers.countPacketBytes) '*uint8=>uint8'], 0);
        if ~isempty(commentIndices)
            [NEV.Data.Comments.TimeStamp, orderOfTS] = sort(Timestamp(commentIndices));
            NEV.Data.Comments.TimeStampSec = double(NEV.Data.Comments.TimeStamp)/double(NEV.MetaTags.TimeRes);
            tempCharSet = tRawData(timeStampBytes+3, commentIndices);
            NEV.Data.Comments.CharSet = tempCharSet(orderOfTS); clear tempCharSet;
            colorFlag = tRawData(timeStampBytes+4, commentIndices);
            NEV.Data.Comments.TimeStampStarted = tRawData(timeStampBytes+5:timeStampBytes+8, commentIndices);
            tempTimeStampStarted = typecast(NEV.Data.Comments.TimeStampStarted(:), 'uint32').';
            NEV.Data.Comments.TimeStampStarted = tempTimeStampStarted(orderOfTS); clear tempTimeStampStarted;
            tempText = char(tRawData(timeStampBytes+9:Trackers.countPacketBytes, commentIndices).');
            NEV.Data.Comments.Text  = tempText(orderOfTS,:); clear tempText;
            
            % Transferring NeuroMotive Events to its own structure
            neuroMotiveEvents = find(NEV.Data.Comments.CharSet == 255);
            NEV.Data.TrackingEvents.TimeStamp = NEV.Data.Comments.TimeStamp(neuroMotiveEvents);
            NEV.Data.TrackingEvents.TimeStampSec = double(NEV.Data.TrackingEvents.TimeStamp)/double(NEV.MetaTags.TimeRes);

            % Parsing NeuroMotive Events
            events = NEV.Data.Comments.Text(neuroMotiveEvents,:);
            for idx = 1:size(events,1)
                splitEvent = strsplit(events(idx,:), ':');
                NEV.Data.TrackingEvents.ROIName{idx} = splitEvent{1};
                NEV.Data.TrackingEvents.ROINum(idx) = str2double(splitEvent{2});
                NEV.Data.TrackingEvents.Event{idx} = splitEvent{3};
                NEV.Data.TrackingEvents.Frame(idx) = str2double(splitEvent{4});
            end
            NEV.Data.Comments.TimeStamp(neuroMotiveEvents) = [];
            NEV.Data.Comments.TimeStampSec(neuroMotiveEvents) = [];
            NEV.Data.Comments.CharSet(neuroMotiveEvents) = [];
            NEV.Data.Comments.TimeStampStarted(neuroMotiveEvents) = [];
            NEV.Data.Comments.TimeStampStartedSec = double(NEV.Data.Comments.TimeStampStarted)/double(NEV.MetaTags.TimeRes);
            NEV.Data.Comments.Text(neuroMotiveEvents,:) = [];
            colorFlag(neuroMotiveEvents) = [];
            
            % Figuring out the text color of the comments that had color
            NEV.Data.Comments.Color = dec2hex(NEV.Data.Comments.TimeStampStarted);
            NEV.Data.Comments.Color(colorFlag == 1,:) = repmat('0', size(NEV.Data.Comments.Color(colorFlag == 1,:)));
            NEV.Data.Comments.TimeStampStarted(colorFlag == 0) = NEV.Data.Comments.TimeStamp(colorFlag == 0);            
            
            clear commentIndices;
        end
        if ~isempty(videoSyncPacketIDIndices)
            NEV.Data.VideoSync.TimeStamp       = Timestamp(videoSyncPacketIDIndices);
            NEV.Data.VideoSync.FileNumber      = tRawData(timeStampBytes+3:timeStampBytes+4, videoSyncPacketIDIndices);
            NEV.Data.VideoSync.FileNumber      = typecast(NEV.Data.VideoSync.FileNumber(:), 'uint16').';
            NEV.Data.VideoSync.FrameNumber     = tRawData(timeStampBytes+5:timeStampBytes+8, videoSyncPacketIDIndices);
            NEV.Data.VideoSync.FrameNumber     = typecast(NEV.Data.VideoSync.FrameNumber(:), 'uint32').';
            NEV.Data.VideoSync.ElapsedTime     = tRawData(timeStampBytes+9:timeStampBytes+12, videoSyncPacketIDIndices);
            NEV.Data.VideoSync.ElapsedTime     = typecast(NEV.Data.VideoSync.ElapsedTime(:), 'uint32').';
            NEV.Data.VideoSync.SourceID        = tRawData(timeStampBytes+13:timeStampBytes+16, videoSyncPacketIDIndices);
            NEV.Data.VideoSync.SourceID        = typecast(NEV.Data.VideoSync.SourceID(:), 'uint32').';
            clear videoSyncPacketIDIndices;
        end
        if ~isempty(trackingPacketIDIndices)
            tmp.TimeStamp     = Timestamp(trackingPacketIDIndices);
            tmp.TimeStampSec  = double(tmp.TimeStamp)/30000;
            % This portion is commented out because it does not contain any
            % information as of yet.
            tmp.ParentID      = tRawData(timeStampBytes+3:timeStampBytes+4, trackingPacketIDIndices);
            tmp.ParentID      = typecast(tmp.ParentID(:), 'uint16').';
            tmp.NodeID        = tRawData(timeStampBytes+5:timeStampBytes+6, trackingPacketIDIndices);
            tmp.NodeID        = typecast(tmp.NodeID(:), 'uint16').';
            tmp.NodeCount     = tRawData(timeStampBytes+7:timeStampBytes+8, trackingPacketIDIndices);
            tmp.NodeCount     = typecast(tmp.NodeCount(:), 'uint16').';
            tmp.MarkerCount   = tRawData(timeStampBytes+9:timeStampBytes+10, trackingPacketIDIndices);
            tmp.MarkerCount   = typecast(tmp.MarkerCount(:), 'uint16').';
            
            tmp.rigidBodyPoints = tRawData(timeStampBytes+11:NEV.MetaTags.PacketBytes, trackingPacketIDIndices);
            tmp.rigidBodyPoints = reshape(typecast(tmp.rigidBodyPoints(:), 'uint16'), size(tmp.rigidBodyPoints, 1)/2, size(tmp.rigidBodyPoints, 2));
            
            if (isfield(NEV, 'ObjTrackInfo'))
                for IDX = 1:size(NEV.ObjTrackInfo,2)
                    emptyChar = find(NEV.ObjTrackInfo(IDX).TrackableName == 0, 1);
                    NEV.ObjTrackInfo(IDX).TrackableName(emptyChar:end) = [];
                    if ~(~isempty(strfind(NEV.ObjTrackInfo(IDX).TrackableName, '1')) || ...
                        ~isempty(strfind(NEV.ObjTrackInfo(IDX).TrackableName, '2')) || ...
                        ~isempty(strfind(NEV.ObjTrackInfo(IDX).TrackableName, '3')) || ...
                        ~isempty(strfind(NEV.ObjTrackInfo(IDX).TrackableName, '4')))
                        nameLength = min(length(NEV.ObjTrackInfo(IDX-1).TrackableName(1:end-1)), length(NEV.ObjTrackInfo(IDX).TrackableName(1:end-1)));
                        if ~strcmpi(NEV.ObjTrackInfo(IDX-1).TrackableName(1:nameLength-1), NEV.ObjTrackInfo(IDX).TrackableName(1:nameLength-1))
                            objectIndex = 1;
                        else
                            objectIndex = objectIndex + 1;
                        end
                        NEV.ObjTrackInfo(IDX).TrackableName(emptyChar) = num2str(objectIndex);
                    end
                    indicesOfEvent = find(tmp.NodeID == IDX-1);
                    if ~isempty(indicesOfEvent)
                        NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).TimeStamp = tmp.TimeStamp(indicesOfEvent);
                        NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).TimeStampSec = tmp.TimeStampSec(indicesOfEvent);
                        NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).ParentID =    tmp.ParentID(indicesOfEvent);
                        NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).NodeCount =   tmp.NodeCount(indicesOfEvent);
                        NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).MarkerCount = tmp.MarkerCount(indicesOfEvent);
                        NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).MarkerCoordinates(size(NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).TimeStamp,2)).X = [];
                        NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).MarkerCoordinates(size(NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).TimeStamp,2)).X = [];
                        for xyIDX = 1:size(NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).TimeStamp,2)
                            NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).MarkerCoordinates(xyIDX).X = ...
                                tmp.rigidBodyPoints(1:2:NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).MarkerCount(xyIDX)*2, indicesOfEvent(xyIDX));
                            NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).MarkerCoordinates(xyIDX).Y = ...
                                tmp.rigidBodyPoints(2:2:NEV.Data.Tracking.(NEV.ObjTrackInfo(IDX).TrackableName).MarkerCount(xyIDX)*2, indicesOfEvent(xyIDX));
                        end
                    end
                end
            end
            clear trackingPacketIDIndices tmp;
        end
        if ~isempty(patientTrigPacketIDIndices)
            NEV.Data.PatientTrigger.TimeStamp    = Timestamp(patientTrigPacketIDIndices);
            NEV.Data.PatientTrigger.TriggerType  = tRawData(timeStampBytes+3:timeStampBytes+4, patientTrigPacketIDIndices);
            NEV.Data.PatientTrigger.TriggerType  = typecast(NEV.Data.PatientTrigger.TriggerType(:), 'uint16').';
            clear patientTrigPacketIDIndices;
        end
        if ~isempty(reconfigPacketIDIndices)
            NEV.Data.Reconfig.TimeStamp     = Timestamp(reconfigPacketIDIndices);
            NEV.Data.Reconfig.ChangeType    = tRawData(timeStampBytes+3:timeStampBytes+4, reconfigPacketIDIndices);
            NEV.Data.Reconfig.ChangeType    = typecast(NEV.Data.Reconfig.ChangeType(:), 'uint16').';
            NEV.Data.Reconfig.CompName      = char(tRawData(timeStampBytes+5:timeStampBytes+20, reconfigPacketIDIndices));
            NEV.Data.Reconfig.ConfigChanged = char(tRawData(timeStampBytes+21:Trackers.countPacketBytes, reconfigPacketIDIndices));
            clear reconfigPacketIDIndices;
        end
        if ~isempty(logEventPacketIDIndices)
            NEV.Data.LogEvent.TimeStamp     = Timestamp(logEventPacketIDIndices);
            tmp.Mode                        = tRawData(timeStampBytes+3:timeStampBytes+4, logEventPacketIDIndices);
            NEV.Data.LogEvent.Mode          = typecast(tmp.Mode(:), 'uint16').';
            NEV.Data.LogEvent.Application   = char(tRawData(timeStampBytes+5:timeStampBytes+20, logEventPacketIDIndices).');
        end
        if ~isempty(recStartPacketIDIndices)
            NEV.Data.RecStartTimes.TimeStamp     = Timestamp(recStartPacketIDIndices);
        end
    end % end if ~isempty(allExtraDataPacketIndices)

    clear Timestamp tRawData count idx;
      
   % now read waveform
    fseek(FID, Trackers.fExtendedHeader + 8, 'bof'); % Seek to location of spikes
    fseek(FID, (Trackers.readPackets(1)-1) * Trackers.countPacketBytes, 'cof');
    NEV.Data.Spikes.WaveformUnit = Flags.waveformUnits;
    NEV.Data.Spikes.Waveform = fread(FID, [(Trackers.countPacketBytes-8)/2 Trackers.readPackets(2)], ...
        [num2str((Trackers.countPacketBytes-8)/2) '*int16=>int16'], 8);
    NEV.Data.Spikes.Waveform(:, [digserIndices allExtraDataPacketIndices]) = []; 

    clear allExtraDataPacketIndices;
    if strcmpi(Flags.waveformUnits, 'uv')
        elecDigiFactors = double(1000./[NEV.ElectrodesInfo(NEV.Data.Spikes.Electrode).DigitalFactor]);
        NEV.Data.Spikes.Waveform = bsxfun(@rdivide, double(NEV.Data.Spikes.Waveform), elecDigiFactors);
        if strcmpi(Flags.WarningStat, 'warning')
            fprintf(1,'\nThe spike waveforms are in unit of uV.\n');
            fprintf(2,'WARNING: This conversion may lead to loss of information.');
            fprintf(1,'\nRefer to help for more information.\n');
        end
    end
end
clear digserIndices;

%% Parse digital data if requested
if ~isempty(DigiValues)
    if strcmpi(Flags.ParseData, 'parse')
        try
            DigiValues = char(DigiValues);
            Inputs                     = {'Digital'; 'AnCh1'; 'AnCh2'; 'AnCh3'; 'AnCh4'; 'AnCh5'; 'PerSamp'; 'Serial'};
            AsteriskIndices   = find(DigiValues == '*');
            DataBegTimestamp = digserTimestamp(AsteriskIndices);
            splitDigiValues = regexp(DigiValues(2:end), '*', 'split')';
            for idx = 1:length(splitDigiValues)
                try
                    if isempty(find(splitDigiValues{idx} == ':', 1))
                        splitDigiValues{idx}(find(splitDigiValues{idx} == '#')) = [];
                        NEV.Data.SerialDigitalIO(idx).Value = splitDigiValues{idx};
                        NEV.Data.SerialDigitalIO(idx).Type = 'Marker';
                    else
                        [tempParsedCommand error] = parseCommand(splitDigiValues{idx});
                        if ~error
                            pcFields = fields(tempParsedCommand);
                            NEV.Data.SerialDigitalIO(idx).Value = splitDigiValues{idx};
                            for fidx = 1:length(pcFields)
                                NEV.Data.SerialDigitalIO(idx).(pcFields{fidx}) = tempParsedCommand.(pcFields{fidx});
                            end
                        else
                            NEV.Data.SerialDigitalIO(idx).Value = splitDigiValues{idx};
                            NEV.Data.SerialDigitalIO(idx).Type = 'UnparsedData';
                            Flags.UnparsedDigitalData = 1;
                        end
                    end
                catch
                    disp(['Error parsing: ' splitDigiValues{idx}]);
                    disp('Please refer to the help for more information on how to properly format the digital data for parsing.');
                end
            end
            % Populate the NEV structure with Timestamp and inputtypes for the
            % digital data
            if ~isempty(DataBegTimestamp)
                c = num2cell(DataBegTimestamp); [NEV.Data.SerialDigitalIO(1:length(NEV.Data.SerialDigitalIO)).TimeStamp] = deal(c{1:end});
                c = num2cell(DataBegTimestamp/NEV.MetaTags.SampleRes); [NEV.Data.SerialDigitalIO.TimeStampSec] = deal(c{1:end});
                c = {Inputs{NEV.Data.SerialDigitalIO.InsertionReason(AsteriskIndices)}}; [NEV.Data.SerialDigitalIO.InputType] = deal(c{1:end});
            end
            clear Inputs DigiValues digserTimestamp;
        catch
            disp(lasterr);
            disp('An error occured during reading digital data. This is due to a problem with formatting digital data.');
            disp('Refer to help ''help openNEV'' for more information on how to properly format the digital data.');
            disp('Try using openNEV with ''noparse'', i.e. openNEV(''noparse'').');
        end
    else
        NEV.Data.SerialDigitalIO.TimeStamp = digserTimestamp;
        NEV.Data.SerialDigitalIO.TimeStampSec = double(digserTimestamp)/30000;
        NEV.Data.SerialDigitalIO.UnparsedData = DigiValues;
        if strcmpi(Flags.Direct, 'direct')
            % Finding the members that have bit 16 as the strobe high
            DShighs = find(NEV.Data.SerialDigitalIO.UnparsedData >= bin2dec('1000000000000000'));
            uniqueDShighs = DShighs([1; find(diff(DShighs)>1)+1]);
            DShighUniqueBin = dec2bin(NEV.Data.SerialDigitalIO.UnparsedData(uniqueDShighs));
            DShighUniqueDec = bin2dec(DShighUniqueBin(:,2:16));
            % Removing the non-strobed-high values from SerialDigitalIO
            extraMembers = setxor(uniqueDShighs, 1:length(NEV.Data.SerialDigitalIO.UnparsedData));  
            NEV.Data.SerialDigitalIO.TimeStamp(extraMembers) = [];
            NEV.Data.SerialDigitalIO.TimeStampSec(extraMembers) = [];
            NEV.Data.SerialDigitalIO.UnparsedData = DShighUniqueDec;
            clear DShighs DShighUniqueBin DShighUniqueDec extraMembers;
        end
        clear digserTimestamp DigiValues
    end
else
    if strcmpi(Flags.ReadData, 'read')
        if strcmpi(Flags.Report, 'report')
            disp('No digital data to read.');
        end
    end
end

if strcmpi(Flags.ParseData, 'parse')
    if Flags.UnparsedDigitalData && strcmpi(Flags.WarningStat, 'warning')
        fprintf(2, 'WARNING: The NEV file contains unparsed digital data.\n');
    end
end

%% Show a report if 'report' is passed as an argument
if strcmpi(Flags.Report, 'report')
    % Displaying report
    disp( '*** FILE INFO **************************');
    disp(['File Name           = ' NEV.MetaTags.Filename]);
    disp(['Filespec            = ' NEV.MetaTags.FileSpec]);
    disp(['Data Duration (min) = ' num2str(round(NEV.MetaTags.DataDuration/NEV.MetaTags.SampleRes/60))]);
    disp(['Packet Counts       = ' num2str(Trackers.countDataPacket)]);
    disp(' ');
    disp( '*** BASIC HEADER ***********************');    
    disp(['Sample Resolution   = ' num2str(NEV.MetaTags.SampleRes)]);
    disp(['Date and Time       = '         NEV.MetaTags.DateTime]);
    disp(['Comment             = '         NEV.MetaTags.Comment(1:64)   ]);
    disp(['                      '         NEV.MetaTags.Comment(65:128) ]);
    disp(['                      '         NEV.MetaTags.Comment(129:192)]);
    disp(['                      '         NEV.MetaTags.Comment(193:256)]);
    disp(['The load time was for NEV file was ' num2str(toc, '%0.1f') ' seconds.']);
end

%% Saving the NEV structure as a MAT file for easy access
if strcmpi(Flags.SaveFile, 'save')
    if exist(matPath, 'file') == 2 && strcmpi(Flags.Overwrite, 'nooverwrite')
        if strcmpi(Flags.WarningStat, 'warning')    
            disp(['File ' matPath ' already exists.']);
            overWrite = input('Would you like to overwrite (Y/N)? ', 's');
        else
            overWrite = 'n';
        end
        if strcmpi(overWrite, 'y')
            if strcmpi(Flags.WarningStat, 'warning')  
                disp('Saving MAT file. This may take a few seconds...');
            end
            save(matPath, 'NEV', '-v7.3');
        else
            if strcmpi(Flags.WarningStat, 'warning')  
                disp('File was not overwritten.');
            end
        end
    elseif exist(matPath, 'file') == 2 && strcmpi(Flags.Overwrite, 'overwrite')
        if strcmpi(Flags.WarningStat, 'warning')
            disp(['File ' matPath ' already exists.']);
            disp('Overwriting the old MAT file. This may take a few seconds...');
        end
        save(matPath, 'NEV', '-v7.3');        
    else
        if strcmpi(Flags.WarningStat, 'warning')
            disp('Saving MAT file. This may take a few seconds...');
        end
        save(matPath, 'NEV', '-v7.3');
    end
    clear overWrite;
end

NEV = killUnwantedChannels(NEV, Flags.selChannels);

if ~nargout
    assignin('base', 'NEV', NEV);
else
    varargout{1} = NEV;
end

fclose(FID);
clear Flags Trackers FID matPath;

function NEV = killUnwantedChannels(NEV, selectedChannels)

if ~strcmpi(selectedChannels, 'all')
    if any(selectedChannels < 1)
        disp('Invalid channel. Channels cannot be 0 or negative values. Channel selection ignored.');
    else
        uniqueChannels = unique(NEV.Data.Spikes.Electrode);
        if ~isempty(setdiff(selectedChannels, uniqueChannels))
            disp('Some of the selected channels in c:xxx command are not in the data file. These will not be loaded.')
        end
        unWantedChannels = setdiff(uniqueChannels, selectedChannels);
        for idx = 1:length(unWantedChannels)
            NEV.Data.Spikes.Waveform(:, NEV.Data.Spikes.Electrode == unWantedChannels(idx)) = [];
            NEV.Data.Spikes.Unit(NEV.Data.Spikes.Electrode == unWantedChannels(idx)) = [];
            NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode == unWantedChannels(idx)) = [];
            NEV.Data.Spikes.Electrode(NEV.Data.Spikes.Electrode == unWantedChannels(idx)) = [];
        end
    end
end

