function varargout = openNSx(varargin)

% openNSx
%
% Opens and reads an NSx file then returns all file information in a NSx
% structure. Works with File Spec 2.1, 2.2, 2.3, and 3.0.
% 
% OUTPUT = openNSx('ver')
% OUTPUT = openNSx(FNAME, 'read', 'report', 'e:xx:xx', 'c:xx:xx', 't:xx:xx', MODE, 'precision', 'skipfactor', 'nozeropad').
% 
% All input arguments are optional. Input arguments can be in any order.
%
%   'ver':        Immediately return version information for openNSx
%                 without processing any files.
%
%   FNAME:        Path of the file to be opened. If FNAME is omitted, a
%                 file selection dialog box will appear.
%
%   'noread':     Do not read the data contained in the file. Return only
%                 header information. ('read' input still accepted for
%                 legacy purposes, but is redundant with default behavior.)
%                 DEFAULT: 'read'
%
%   'report':     Show a summary report if user passes this argument.
%                 DEFAULT: No report.
%
%   'electrodes',XX
%   'e:XX':       User can specify which electrodes need to be read. Values 
%                 provided for XX represent one or more electrode IDs,
%                 which are only available in the *.cmp file (map file)
%                 accompanying the Blackrock electrode array product. These
%                 are NOT the same as the channel IDs stored in
%                 `NSx.ElectrodesInfo.ChannelID`; see the file specification
%                 documentation for more information. XX can take the form
%                 of any string interpretable by MATLAB `num2str`, e.g., 
%                 '1:10' or '[1 3 5]', or for the separate key-value form,
%                 the actual numerical values. The requested electrode IDs
%                 must lie in the set of electrode IDs corresponding to 
%                 channels of data in the data file. Note that when the set
%                 of channels requested is not continguous, `openNSx` still
%                 reads the contiguous block of channels (min:max), then
%                 subselects afterward to produce the requested output. 
%                 This may result in higher peak memory usage than 
%                 expected. Rows of data in the output retain ordering in
%                 XX, e.g. for XX=[X1 X2 X3 ...], the rows of NSx.Data will
%                 correspond to X1, X2, X3, etc., even if X1, X2, X3 are
%                 not sorted in ascending order. Use of this option
%                 requires that the user also provide the *.cmp mapfile
%                 provided by Blackrock (when prompted), and that 
%                 KTUEAMapFile is present in path.
%                 DEFAULT: will read all existing electrodes.
%
%   'channels',XX
%   'c:XX':       User can specify which channels need to be read. Values 
%                 provided for XX represent one or more indices into the
%                 set of channels contained in the data file. These
%                 are NOT the same as the channel IDs stored in
%                 `NSx.ElectrodesInfo.ChannelID`; see the file specification
%                 documentation for more information. XX can take the form
%                 of any string interpretable by MATLAB `num2str`, e.g., 
%                 '1:10' or '[1 3 5]', or for the separate key-value form,
%                 the actual numerical values. The requested channel
%                 indices must lie in the range from 1 to the number of
%                 channels in the data file. Note that when the set of
%                 channels requested is not continguous, `openNSx` still
%                 reads the contiguous block of channels (min:max), then
%                 subselects afterward to produce the requested output. 
%                 This may result in higher peak memory usage than 
%                 expected. Channels in the output retain ordering in XX,
%                 e.g. for XX=[X1 X2 X3 ...], the rows of NSx.Data will
%                 correspond to X1, X2, X3, etc., even if X1, X2, X3 are
%                 not sorted in ascending order.
%                 DEFAULT: will read all existing analog channels.
%
%   'duration',XX
%   't:XX':       Specify the beginning and end of the data to be read. 
%                 Values in XX are interpreted as though all data samples
%                 in the file were recorded contiguously from time 0,
%                 without regard for actual segment timestamps or durations.
%                 The units of XX are determined by the MODE input below,
%                 and converted into units of samples of data. Thus, while
%                 duration may be specified with units of time, such as
%                 seconds, XX will not be interpreted as "real" time. See
%                 below for an example of the implication of this behavior.
%                 However, the data and metadata returned by `openNSx` will
%                 reflect the correct segment boundaries and timestamps.
%                 Note that timestamps are provided in units of
%                 `NSx.MetaTags.TimeRes`, not the sampling rate of the
%                 file. XX can take the form of any string interpretable
%                 by MATLAB `num2str`, e.g., '1:10' or '[1 3 5]', or for
%                 the separate key-value form, the actual numerical values.
%                 If the start time is before the start of the file, it
%                 will be changed to the start of the file with a warning.
%                 If it is past the end of the file, `openNSx` will exit
%                 with an error message. If the end time is greater than
%                 the length of data the user will be prompted to use the
%                 last datapoint in the file with a warning. 
%                 DEFAULT: will read the entire file.
%                 EXAMPLE: In two-NSP recordings with a clock restart for
%                 syncing, there is typically a short data segment, with
%                 samples captured before the clock restart, followed by
%                 a longer segment. The duration of the first segment
%                 is usually longer than the timestamp of the second
%                 segment would suggest, e.g., segment timestamps of
%                 [0 116] but durations of [8000 100000] (for sampling rate
%                 of 2,0000 samples/sec and TimestampTimeResolution of
%                 30000). The values in XX will be interpreted as though
%                 the file contained a single data segment with timestamp
%                 0 and duration 108000. However, if the requested data
%                 falls across the real segment boundaries, the data will
%                 be appropriately grouped in cell arrays, and the MetaTags
%                 timestamps will reflect the appropriate values.
%
%   MODE:         Specify the units of duration values specified with 
%                 'duration' (or 't:XX:YY') input. Valid values of MODE are
%                       'sec', 'secs', 'second', 'seconds'
%                       'min', 'mins', 'minute', 'minutes'
%                       'hour', 'hours'
%                       'sample', 'samples'
%                 Note that when MODE is 'samples', duration input (see
%                 above) must be greater than or equal to 1; however, for
%                 all other values of MODE, the origin is 0.
%                 DEFAULT: 'sample'
%
%   'uV':         Read the recording waveforms in unit of uV instead of raw
%                 values. Note that this conversion requires 'double'
%                 precision; if this argument is provided and precision has
%                 been set to 'int16' or 'short', it will be updated to 
%                 'double' with a warning.
%
%   'precision',P
%   'p:P':
%   P             Specify the precision P for data read from the NSx file.
%                 Valid options are 'double' (64-bit floating point) or
%                 'int16' (or, equivalently, 'short'; 16-bit signed
%                 integer). Data are stored in the file as int16 values. 
%                 Still, while 'int16' uses less memory, be mindful of
%                 the limitations of integer data types
%                 (https://www.mathworks.com/help/matlab/numeric-types.html).
%                 Note that if the argument 'uV' is provided (conversion
%                 from raw values to uV units), the precision will be
%                 automatically set to 'double' if it is not already.
%                 DEFAULT: 'int16'.
%
%   'skipfactor',S
%   's:S':        Decimate data read from disk, e.g., to quickly preview
%                 data. The integer S will determine how many samples to
%                 skip. For example, if S is 2 then every other sample is
%                 read. This action is decimation only: no anti-aliasing
%                 filter is applied.
%                 DEFAULT: 1 (every sample read)
%
%   'zeropad':    Prepend the data with zeros to compensate for non-zero 
%                 start time. Note that timestamps in newer data files may
%                 be in the 10^18 range. Prepending this many zeros is not
%                 advisable for normal computer systems, nor does it carry
%                 the same logic as older data files where time zero was
%                 relevant to the recording.
%                 DEFAULT: No zero padding.
%
%   'noalign':    Do not apply bug fix for clock drift in Central release
%                 7.6.0. Only applies to files with precision time protocol
%                 (PTP) nanosecond-resolution timestamps. By default,
%                 samples may be added (by duplication) or removed (by
%                 deletion) to restore clock alignment. Changes are made at
%                 evenly spaced points throughout the data segments within
%                 the file. With alignment enabled, the starting timestamp
%                 of the segment is sufficient to infer the timestamps of
%                 all subsequent samples in the segment. Without alignment,
%                 the timestamps of each sample must be provided to know
%                 their timing (returned in NSx.Time). Note the increase in
%                 memory footprint required for this change.
%                 DEFAULT: Alignment occurs with warnings.
%
%   'nosegment':  Do not segment the file around longer-than-expected
%                 intervals between timestamps (i.e., pauses). Only applies
%                 to files with precision time protocol (PTP)
%                 nanosecond-resolution timestamps. See the documentation
%                 for 'max_tick_multiple' for a more detailed explanation. 
%                 DEFAULT: segment the data around pauses.
%
%   'max_tick_multiple', M:
%                 Newer data files use PTP (precision time protocol) and
%                 timestamp each sample of data, instead of only the first
%                 sample in a frame of contiguous samples. To detect pauses
%                 in a PTP recording, openNSx processes the file in frames:
%                 it reads the timestamp of the first and last packets in
%                 each frame (see 'packets_per_frame') and checks whether
%                 the elapsed time is greater than it should be, assuming
%                 contiguously recorded packets at the expected sampling
%                 rate. The threshold M for the difference of elapsed time
%                 is set as a multiple of the expected sampling interval.
%                 If M is too small, openNSx will detect spurious pauses.
%                 If it is too high, pauses will be missed. Note that due
%                 to jitter in sample timing, this value should be set in
%                 coordination with the number of packets in each frame
%                 (see 'packets_per_frame') to ensure the sum of jittered
%                 sampling intervals does not exceed the detection
%                 threshold.
%                 DEFAULT: 2 (equivalent to missing one sample)
%
%   'packets_per_frame', P:
%                 Newer data files use PTP (precision time protocol) and
%                 timestamp each sample of data, instead of only the first
%                 sample in a frame of contiguous samples. To detect pauses
%                 in a PTP recording, openNSx processes the file in frames,
%                 each containing P packets: it reads the timestamp of the
%                 first and last packets in each frame and checks whether
%                 the elapsed time is greater than it should be, assuming
%                 contiguously recorded packets at the expected sampling
%                 rate (see 'max_tick_multiple'). The number of frames F is
%                 given by CEIL(TOTAL_PACKETS/P), where TOTAL_PACKETS is
%                 the number of packets in the file. Note that this method
%                 reads only F+1 timestamps from disk if there are no
%                 pauses detected. For each frame containing one or more
%                 detected pauses, all P timestamps in the frame are read
%                 from disk to identify the specific samples between which
%                 the pauses occur. Thus, P can be increased to lower F,
%                 but it should not be so large that a vector of P doubles
%                 would not fit in memory. Note also that because of jitter
%                 in sample timing, setting this value too large may lead
%                 to spurious detections (i.e., the sum of jitter could be
%                 greater than the detection threshold).
%                 DEFAULT: 100,000 packets per frame.
%
%   OUTPUT:       The NSx structure.
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
%   version of the datafile will be read, where only every 5th sample is
%   read.
%
%   Example 2:
%   openNSx('read','c:15:30');
%
%   In the example above, the user will be prompted for the file. The file
%   will be read using 'int16' precision as default. All time points of
%   channels 15 through 30 will be read.
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
% 7.1.0.0: April 14, 2020
%   - Added option to load the data without zero padding to compensate for
%     a non-zero start time. (David Kluger)
%   - Bug fixes and documentation updates (David Kluger)
%
% 7.1.1.0: June 11, 2020
%   - Fixed a bug related to fread and MATLAB 2020a.
%
% 7.3.0.0: September 11, 2020
%   - Fixed a bug related to fread and MATLAB 2020a.
%   - Gives a warning about FileSpec 3.0 and gives the user options for how
%     to proceed.
%   - Added a warning about the data unit and that by default it in the
%     unit of 250 nV or 1/4 µV.
%   - If the units are in "raw", ths correct information is now written to
%     the electrodes header: 250 nV (raw).
%
% 7.3.1.0: October 2, 2020
%   - If the units are in µV (openNSx('uv'), ths correct information is now
%     written to the electrodes header: 1000 nV (raw).
%
% 7.3.2.0: October 23, 2020
%   - Fixed a typo.
%
% 7.4.0.0: October 29, 2020
%   - Undid changes made to AnalogUnit and instead implemented
%     NSx.ElectrodesInfo.Resolution to show what the resolution of the data
%     is. By default, the resolution is set to 0.250 µV. If used with
%     parameter 'uv', the resolution will be 1 µV. To always convert the
%     data to µV, divide NSx.Data(CHANNEL,:) by
%     NSx.ElectrodesInfo(CHANNEL).Resolution.
%
% 7.4.1.0: April 20, 2021
%   - Fixed a bug related to file opening.
%
% 7.4.2.0: May 5, 2021
%   - Fixed a bug related to NeuralSG file format (File Spec 2.1).
%
% 7.4.3.0: July 16, 2021
%   - Fixed a minor bug for when the data header is not written properly
%     and the data needs to be used to calculate the data length.
%
% 7.4.4.0: April 1, 2023
%   - Accounts for many segments in files for clock drift correction
%   - Changed 'zeropad' default behavior to be 'no'
%
% 7.4.5.0: October 5, 2023
%   - Bank numbers on new files are not alpha which causes problems on save
%
% 7.4.6.0: December 6, 2023
%   - Better support for reading files recorded from Gemini systems
%   - Improved speed and memory usage for Gemini system recordings
%   - Change messages about errors to actual errors
%   - NPMK SettingsManager, getFile, and NPMKverChecker made optional
%   - Force 'double' precision (with warning) if conversion to uV requested
%   - Repair skipfactor implementation
%   - Clean up documentation
%   - Clean up code
%
% 7.4.6.1: January 30, 2024
%   - Bug fix: mishandling of numerical input arguments
%   - Bug fix: support noncontiguous channel output
%
% 7.4.6.2: April 26, 2024
%   - Add feature to disable data segmentation for file spec >=3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define the NSx data structure and sub-branches.
NSx          = struct('MetaTags',[],'Data',[],'Time',[],'RawData', []);
NSx.MetaTags = struct('FileTypeID',[],'SamplingLabel',[],'ChannelCount',[],'SamplingFreq',[], 'TimeRes', [], ...
                      'ChannelID',[],'DateTime',[],'DateTimeRaw',[], 'Comment', [], 'FileSpec', [], ...
                      'Timestamp', [], 'DataPoints', [], 'DataDurationSec', [], 'openNSxver', [], 'Filename', [], 'FilePath', [], ...
                      'FileExt', []);

NSx.MetaTags.openNSxver = '7.4.6.2';

%% Check for the latest version of NPMK
if exist('NPMKverChecker','file')==2
    NPMKverChecker
end

%% Define constants and defaults
extHeaderEntrySize = 66;
NSx.RawData.PausedFile = 0;
syncShift = 0;
flagFoundSettingsManager = exist('settingsManager','file')==2;
flagFoundGetFile = exist('getFile','file')==2;
NPMKSettings = [];
if flagFoundSettingsManager
    NPMKSettings = settingsManager;
end

% Default values
flagReport = 0;
flagReadData = 1;
flagModifiedTime = 0;
flagMultiNSP = 1;
flagZeroPad = 0;
flagAlign = 1;
flagSegment = 1;
flagConvertToUv = 0;
flagOneSamplePerPacket = 0;
requestedTimeScale = 'sample';
requestedPrecisionType = 'int16';
requestedSkipFactor = 1;
requestedPacketsPerFrame = 100000;
requestedMaxTickMultiple = 2;
requestedChannelIDs = [];
requestedChannelIndex = [];
requestedFileName = '';

%% Process input arguments
next = '';
for i=1:length(varargin)
    inputArgument = varargin{i};
    if ischar(inputArgument) && strcmpi(inputArgument, 'ver')
        varargout{1} = NSx.MetaTags.openNSxver;
        return;
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'channels')
        next = 'channels';
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'skipfactor')
        next = 'skipfactor';
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'electrodes')
        next = 'electrodes';
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'duration')
        next = 'duration';
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'precision')
        next = 'precision';
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'packets_per_frame')
        next = 'packets_per_frame';
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'max_tick_multiple')
        next = 'max_tick_multiple';
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'report')
        flagReport = 1;
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'noread')
        flagReadData = 0;
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'nomultinsp')
        flagMultiNSP = 0;
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'zeropad')
        flagZeroPad = 1;
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'uV')
        flagConvertToUv = 1;
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'noalign')
        flagAlign = false;
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'nosegment')
        flagSegment = false;
    elseif ischar(inputArgument) && strcmpi(inputArgument, 'read')
        flagReadData = 1;
    elseif (ischar(inputArgument) && strncmp(inputArgument, 't:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'duration')
        if strncmp(inputArgument, 't:', 2)
            inputArgument(1:2) = [];
            inputArgument = str2num(inputArgument); %#ok<ST2NM>
        elseif ischar(inputArgument)
            inputArgument = str2num(inputArgument); %#ok<ST2NM>
        end
        assert(isnumeric(inputArgument),'Must provide duration input that can be resolved to a numeric value');
        flagModifiedTime = 1;
        requestedStartValue = inputArgument(1);
        requestedEndValue = inputArgument(end);
        next = '';
    elseif (ischar(inputArgument) && strncmp(inputArgument, 'e:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'electrodes')
        assert(exist('KTUEAMapFile','file')==2,'To read data by ''electrodes'' the function KTUEAMapFile needs to be in path.');
        mapFile = KTUEAMapFile;
        if strncmp(inputArgument, 'e:', 2)
            requestedElectrodes = str2num(inputArgument(3:end)); %#ok<ST2NM>
        elseif ischar(inputArgument)
            requestedElectrodes = str2num(inputArgument); %#ok<ST2NM>
        else
            requestedElectrodes = inputArgument;
        end
        assert(isnumeric(inputArgument),'Must provide electrode input that can be resolved to a numeric value');
        if min(requestedElectrodes)<1 || max(requestedElectrodes)>128
            assert(min(requestedElectrodes)>=1 && max(requestedElectrodes)<=128, 'The electrode number cannot be less than 1 or greater than 128.');
        end
        requestedChannelIDs = nan(1,length(requestedElectrodes));
        for chanIDX = 1:length(requestedElectrodes)
            requestedChannelIDs(chanIDX) = mapFile.Electrode2Channel(requestedElectrodes(chanIDX));
        end
        next = '';
    elseif (ischar(inputArgument) && strncmp(inputArgument, 's:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'skipFactor')
        if strncmp(inputArgument, 's:', 2)
            requestedSkipFactor = str2num(inputArgument(3:end)); %#ok<ST2NM>
        elseif ischar(inputArgument)
            requestedSkipFactor = str2num(inputArgument); %#ok<ST2NM>
        else
            requestedSkipFactor = inputArgument;
        end
        assert(isnumeric(requestedSkipFactor) && isscalar(requestedSkipFactor) && floor(requestedSkipFactor)==requestedSkipFactor,'Must provide skipfactor input that can be resolved to a scalar integer value')
        next = '';
    elseif (ischar(inputArgument) && strncmp(inputArgument, 'c:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/') || strcmpi(next, 'channels')
        if strncmp(inputArgument, 'c:', 2)
            requestedChannelIndex = str2num(inputArgument(3:end)); %#ok<ST2NM>
        elseif ischar(inputArgument)
            requestedChannelIndex = str2num(inputArgument); %#ok<ST2NM>
        else
            requestedChannelIndex = inputArgument;
        end
        assert(isnumeric(requestedChannelIndex),'Must provide channel input that can be resolved to a numeric value');
        next = '';
    elseif (ischar(inputArgument) && any(strcmpi(inputArgument,{'double','int16','short'})) || (strncmp(varargin{i}, 'p:', 2) && inputArgument(3) ~= '\' && inputArgument(3) ~= '/')) || strcmpi(next, 'precision')
        if strncmpi(inputArgument, 'p:', 2)
            precisionTypeRaw = inputArgument(3:end);
        else
            precisionTypeRaw = inputArgument;
        end
        switch precisionTypeRaw
            case {'int16','short'}
                requestedPrecisionType = 'int16';
            case 'double'
                requestedPrecisionType = 'double';
            otherwise
                error('Precision type is not valid. Refer to ''help'' for more information.');
        end
        next = '';
    elseif strcmpi(next, 'packets_per_frame')
        if ischar(inputArgument)
            requestedPacketsPerFrame = str2double(inputArgument);
        else
            requestedPacketsPerFrame = inputArgument;
        end
    elseif strcmpi(next, 'max_tick_multiple')
        if ischar(inputArgument)
            requestedMaxTickMultiple = str2double(inputArgument);
        else
            requestedMaxTickMultiple = inputArgument;
        end
    elseif ischar(inputArgument) && ...
            (strncmpi(inputArgument,'hours',4) || strncmpi(inputArgument,'hrs',2) || ...
            strncmpi(inputArgument,'minutes',3) || strncmpi(inputArgument,'mins',3) || ...
            strncmpi(inputArgument,'seconds',3) || strncmpi(inputArgument,'secs',3) || ...
            strncmpi(inputArgument,'samples',4))
        requestedTimeScale = inputArgument;
    elseif ischar(inputArgument) && length(inputArgument)>3 && ...
            (strcmpi(inputArgument(3),'\') || ...
            strcmpi(inputArgument(1),'/') || ...
            strcmpi(inputArgument(2),'/') || ...
            strcmpi(inputArgument(1:2), '\\') || ...
            strcmpi(inputArgument(end-3), '.'))
        requestedFileName = inputArgument;
        assert(exist(requestedFileName, 'file')==2,'The file does not exist.');
    else
        error(['Invalid argument ''' inputArgument '''.']);
    end
end
clear next;

% check uV conversion versus data type
if flagReadData && flagConvertToUv && ~strcmpi(requestedPrecisionType,'double')
    warning("Conversion to uV requires double precision; overriding user request '%s' to comply",requestedPrecisionType);
    requestedPrecisionType = 'double';
end

% warn if only reading header information
if ~flagReadData
    warning('Reading the header information only.');
end

% start the report if requested
if flagReport
    disp(['openNSx ' NSx.MetaTags.openNSxver]);
end

%% Identify data file name, path, and extension
%  for later use, and validate the entry.
if isempty(requestedFileName)
    title = 'Choose an NSx file...';
    filterSpec = '*.ns*';
    if flagFoundGetFile
        [requestedFileName, requestedFilePath] = getFile(filterSpec, title);
    else
        [requestedFileName, requestedFilePath] = uigetfile(filterSpec, title);
    end
    assert(ischar(requestedFileName),'No file selected');
    [~, ~, requestedFileExtension] = fileparts(requestedFileName);
else
    if isempty(fileparts(requestedFileName))
        requestedFileName = which(requestedFileName);
    end
    [requestedFilePath,requestedFileName, requestedFileExtension] = fileparts(requestedFileName);
    requestedFileName = [requestedFileName requestedFileExtension];
    requestedFilePath  = [requestedFilePath '/'];
end
assert(ischar(requestedFileName)||requestedFileName~=0,'Could not identify file to read');
fileFullPath = fullfile(requestedFilePath, requestedFileName);
[NSx.MetaTags.FilePath, NSx.MetaTags.Filename, NSx.MetaTags.FileExt] = fileparts(fileFullPath);

% Check to see if 512 setup and calculate offset
if flagMultiNSP
    flag512 = regexp(requestedFileName, '-i[0123]-', 'ONCE');
    if ~isempty(flag512)
        syncShift = multiNSPSync(fullfile(requestedFilePath, requestedFileName));
    else
        flagMultiNSP = 0;
    end
end

%% Loading .x files for multiNSP configuration
if strcmpi(requestedFileExtension(2:4), 'ns6') && length(requestedFileExtension) == 5
    requestedFilePath(1) = requestedFileName(end);
    requestedFileName(end) = [];
end

%% Measure time required to load data
tic;

%% Process file
FID = fopen([requestedFilePath requestedFileName], 'r', 'ieee-le');
try
    
    %% Read Headers
    NSx.MetaTags.FileTypeID = fread(FID, [1,8], 'uint8=>char');
    if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALSG')
        
        %% Read Basic Header
        basicHeaderBytes           = fread(FID, 24, '*uint8');
        NSx.MetaTags.FileSpec      = '2.1';
        NSx.MetaTags.SamplingLabel = char(basicHeaderBytes(1:16));
        NSx.MetaTags.TimeRes       = double(30000);
        NSx.MetaTags.SamplingFreq  = NSx.MetaTags.TimeRes / double(typecast(basicHeaderBytes(17:20),'uint32'));
        channelCount               = double(typecast(basicHeaderBytes(21:24),'uint32'));
        NSx.MetaTags.ChannelCount  = channelCount;
        
        %% Read Extended Header
        extendedHeaderSize = channelCount*4;
        extendedHeaderBytes = fread(FID, extendedHeaderSize, '*uint8');
        NSx.MetaTags.ChannelID = typecast(extendedHeaderBytes, 'uint32');
        try
            t                          = dir(fileFullPath);
            NSx.MetaTags.DateTime      = t.date;
        catch ME2
            warning('openNSx:NEURALSG_datetime','Could not compute date from file: %s',ME2.message)
            NSx.MetaTags.DateTime  = '';
        end
        timestampSize             = 4;
    elseif or(strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD'), strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP'))
        
        %% Read Basic Header
        basicHeaderBytes           = fread(FID, 306, '*uint8');
        NSx.MetaTags.FileSpec      = [num2str(double(basicHeaderBytes(1))) '.' num2str(double(basicHeaderBytes(2)))];
        %BasicHeaderSize            = double(typecast(BasicHeader(3:6), 'uint32'));
        NSx.MetaTags.SamplingLabel = char(basicHeaderBytes(7:22))';
        NSx.MetaTags.Comment       = char(basicHeaderBytes(23:278))';
        NSx.MetaTags.TimeRes       = double(typecast(basicHeaderBytes(283:286), 'uint32'));
        NSx.MetaTags.SamplingFreq  = double(30000 / double(typecast(basicHeaderBytes(279:282), 'uint32')));
        t                          = double(typecast(basicHeaderBytes(287:302), 'uint16'));
        channelCount               = double(typecast(basicHeaderBytes(303:306), 'uint32'));
        NSx.MetaTags.ChannelCount  = channelCount;
        if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD')
            timestampSize = 4;
            timestampType = 'uint32';
        elseif strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP')
            timestampSize = 8;
            timestampType = 'uint64';
        end
        
        % Removing extra garbage characters from the Comment field.
        NSx.MetaTags.Comment(find(NSx.MetaTags.Comment==0,1):end) = 0;
        
        %% Read Extended Header
        extendedHeaderSize = double(channelCount * extHeaderEntrySize);
        extendedHeaderBytes = fread(FID, extendedHeaderSize, '*uint8');
        for headerIDX = 1:channelCount
            byteOffset = double((headerIDX-1)*extHeaderEntrySize);
            NSx.ElectrodesInfo(headerIDX).Type = char(extendedHeaderBytes((1:2)+byteOffset))';
            assert(strcmpi(NSx.ElectrodesInfo(headerIDX).Type, 'CC'),'extended header not supported');
            
            NSx.ElectrodesInfo(headerIDX).ChannelID = typecast(extendedHeaderBytes((3:4)+byteOffset), 'uint16');
            NSx.ElectrodesInfo(headerIDX).Label = char(extendedHeaderBytes((5:20)+byteOffset))';
            NSx.ElectrodesInfo(headerIDX).ConnectorBank = extendedHeaderBytes(21+byteOffset);
            NSx.ElectrodesInfo(headerIDX).ConnectorPin   = extendedHeaderBytes(22+byteOffset);
            NSx.ElectrodesInfo(headerIDX).MinDigiValue   = typecast(extendedHeaderBytes((23:24)+byteOffset), 'int16');
            NSx.ElectrodesInfo(headerIDX).MaxDigiValue   = typecast(extendedHeaderBytes((25:26)+byteOffset), 'int16');
            NSx.ElectrodesInfo(headerIDX).MinAnalogValue = typecast(extendedHeaderBytes((27:28)+byteOffset), 'int16');
            NSx.ElectrodesInfo(headerIDX).MaxAnalogValue = typecast(extendedHeaderBytes((29:30)+byteOffset), 'int16');
            NSx.ElectrodesInfo(headerIDX).AnalogUnits    = char(extendedHeaderBytes((31:46)+byteOffset))';
            if flagConvertToUv
                NSx.ElectrodesInfo(headerIDX).Resolution = 1;
            else
                NSx.ElectrodesInfo(headerIDX).Resolution = ...
                    round(double(NSx.ElectrodesInfo(headerIDX).MaxAnalogValue) / double(NSx.ElectrodesInfo(headerIDX).MaxDigiValue),4);
            end
            %         if strcmpi(waveformUnits, 'uV')
            %             NSx.ElectrodesInfo(headerIDX).AnalogUnits    = '1000 nV (raw)   ';
            %         else
            %             conversion = int16(double(NSx.ElectrodesInfo(headerIDX).MaxAnalogValue) / double(NSx.ElectrodesInfo(headerIDX).MaxDigiValue)*1000);
            %             NSx.ElectrodesInfo(headerIDX).AnalogUnits    = [num2str(conversion), ' nV (raw)    '];
            %         end
            NSx.ElectrodesInfo(headerIDX).HighFreqCorner = typecast(extendedHeaderBytes((47:50)+byteOffset), 'uint32');
            NSx.ElectrodesInfo(headerIDX).HighFreqOrder  = typecast(extendedHeaderBytes((51:54)+byteOffset), 'uint32');
            NSx.ElectrodesInfo(headerIDX).HighFilterType = typecast(extendedHeaderBytes((55:56)+byteOffset), 'uint16');
            NSx.ElectrodesInfo(headerIDX).LowFreqCorner  = typecast(extendedHeaderBytes((57:60)+byteOffset), 'uint32');
            NSx.ElectrodesInfo(headerIDX).LowFreqOrder   = typecast(extendedHeaderBytes((61:64)+byteOffset), 'uint32');
            NSx.ElectrodesInfo(headerIDX).LowFilterType  = typecast(extendedHeaderBytes((65:66)+byteOffset), 'uint16');
        end
        
        % Parse DateTime
        NSx.MetaTags.DateTimeRaw = t.';
        NSx.MetaTags.DateTime = char(datetime(t(1), t(2), t(4), t(5), t(6), t(7)));
    else
        error('Unsupported file spec %s', NSx.MetaTags.FileSpec);
    end
    
    % Check zeropad if timeres is 1e9
    if NSx.MetaTags.TimeRes > 1e5
        assert(~flagZeroPad,'No zeropad for nanosecond resolution timestamps');
    end
    
    % Copy ChannelID to MetaTags for filespec 2.2, 2.3, and 3.0 for compatibility with filespec 2.1
    if or(strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD'), strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP'))
        NSx.MetaTags.ChannelID = [NSx.ElectrodesInfo.ChannelID]';
    end
    
    %% Identify key points in file
    f.EOexH = double(ftell(FID));
    fseek(FID, 0, 'eof');
    f.EOF = double(ftell(FID));
    
    %% Save raw headers for saveNSx
    NSx.RawData.Headers = [uint8(NSx.MetaTags.FileTypeID(:)); basicHeaderBytes(:); extendedHeaderBytes(:)];
    
    %% Verify channel ID and index variables
    if isempty(requestedChannelIndex)
        if isempty(requestedChannelIDs)
            requestedChannelIndex = 1:length(NSx.MetaTags.ChannelID);
        else
            requestedChannelIndex = nan(1,length(requestedChannelIDs));
            for idx = 1:length(requestedChannelIDs)
                assert(ismember(requestedChannelIDs(idx), NSx.MetaTags.ChannelID),'Channel ID %d does not exist in this file',requestedChannelIDs(idx));
                requestedChannelIndex(idx) = find(NSx.MetaTags.ChannelID == requestedChannelIDs(idx),1);
            end
        end
    end
    if isempty(requestedChannelIDs)
        requestedChannelIDs = NSx.MetaTags.ChannelID(requestedChannelIndex);
    end
    assert(all(requestedChannelIndex<=channelCount),'Channel indices must be less than or equal to the total number of channels in the file (%d)',channelCount);
    assert(all(ismember(requestedChannelIDs,NSx.MetaTags.ChannelID)),'Channel IDs must be present in the set of channel IDs in the file');
    NSx.MetaTags.ChannelCount = length(requestedChannelIDs);
    requestedFirstChannel = min(requestedChannelIndex);
    requestedLastChannel = max(requestedChannelIndex);
    numChannelsToRead = length(requestedFirstChannel:requestedLastChannel);
    
    %% Central v7.6.0 needs corrections for PTP clock drift - DK 20230303
    if NSx.MetaTags.TimeRes > 1e5
        packetSize = 1 + timestampSize + 4 + channelCount*2; % byte (Header) + uint64 (Timestamp) + uint32 (Samples, always 1) + int16*nChan (Data)
        numPacketsTotal = floor((f.EOF - f.EOexH)/packetSize);
        fseek(FID, f.EOexH + 1 + timestampSize, 'bof'); % byte (Header) + uint64 (Timestamp)
        patchCheck = fread(FID,10,'uint32',packetSize-4); % read "samples" counts from 10 packets
        if sum(patchCheck) == length(patchCheck) % verify all 1
            flagOneSamplePerPacket = true;
        end
    end

    % either or both align and segment must be enabled
    if ~flagSegment && flagOneSamplePerPacket
        if flagAlign
            flagAlign = false;
            warning('Disabling alignment because segmentation was disabled')
        end
    end
    
    %% Identify and describe data segments
    fseek(FID, f.EOexH, 'bof');
    if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALSG')
        NSx.MetaTags.Timestamp = 0; % No timestamp otherwise
        NSx.MetaTags.DataPoints = double(f.EOF-f.EOexH)/(channelCount*2);
        NSx.MetaTags.DataDurationSec = NSx.MetaTags.DataPoints/NSx.MetaTags.SamplingFreq;
    elseif or(strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD'), strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP'))
        if flagOneSamplePerPacket
            if flagSegment
                
                % Clock drift patch kills ability to segment files. This check will
                % allow segments to be reintroduced into the data structures if a
                % timestamp difference of 200% greater than expected is identified
                fseek(FID,f.EOexH + 1,'bof'); % + byte (header)

                % Process file in frames. initialize with the first packet's
                % timestamp.
                % For each frame, read the timestamp of the last packet.
                % if the difference from the previous frame's last timestamp is
                % larger than expected given consistent sampling rates, define a
                % segment.
                % Move to the next frame.
                ticksPerSample = NSx.MetaTags.TimeRes/NSx.MetaTags.SamplingFreq;
                minimumPauseLength = requestedMaxTickMultiple*ticksPerSample;
                timestampFirst = fread(FID,1,timestampType);
                numPacketsProcessed = 0;
                segmentTimestamps = nan(1,1e3);
                segmentTimestamps(1) = timestampFirst;
                segmentDatapoints = nan(1,1e3);
                segmentDurations = nan(1,1e3);
                currSegment = 1;
                while double(ftell(FID)) < (f.EOF-(packetSize-1-8))
                    
                    % frames have 'packets_per_frame' packets until the end of the
                    % file, when the frame may have fewer packets
                    % number of packets per frame includes first/last packet, which
                    % means there is one fewer gap than the number of packets
                    currPacketStartByte = double(ftell(FID)) - 8 - 1;
                    frameNumPackets = min(requestedPacketsPerFrame, (f.EOF - currPacketStartByte)/packetSize);
                    if abs(round(frameNumPackets)-frameNumPackets)>0.1
                        warning('File not packet-aligned')
                    end
                    bytesToFrameLastTimestamp = packetSize*(frameNumPackets-1) - 8;
                    
                    % compute the ticks expected to elapse in this frame with the
                    % smallest detectable pause (2x sample time, or 66.6 usec)
                    expectedTicksElapsedNoPause = (frameNumPackets-1) * ticksPerSample;
                    expectedTicksElapsedMinPause = expectedTicksElapsedNoPause + (minimumPauseLength - ticksPerSample);
                    
                    % seek to last packet of this frame and read timestamp
                    fseek(FID, bytesToFrameLastTimestamp, 'cof');
                    timestampLast = fread(FID,1,timestampType);
                    
                    % check whether elapsed time for this frame meets or exceeds
                    % expected length with minimum gap
                    actualTicksElapsed = timestampLast - timestampFirst;
                    if actualTicksElapsed >= expectedTicksElapsedMinPause
                        
                        % a gap exists in this frame; we need to identify where it
                        % occurs
                        % save file pointer position
                        currBytePosition = ftell(FID);
                        
                        % rewind to prior last_timestamp
                        fseek(FID, -(bytesToFrameLastTimestamp+8+8), 'cof');
                        
                        % read all timestamps in this frame
                        timestamps = fread(FID, frameNumPackets, timestampType, packetSize-8)';
                        
                        % find gaps and store if found
                        tsDiffs = diff(timestamps);
                        vals = find(tsDiffs > minimumPauseLength);
                        for jj=1:length(vals)
                            numDatapointsLastSegment = numPacketsProcessed - sum(segmentDatapoints(~isnan(segmentDatapoints))) + vals(jj);
                            segmentDatapoints(currSegment) = numDatapointsLastSegment;
                            segmentDurations(currSegment) = timestamps(vals(jj)) - segmentTimestamps(currSegment) + 1;
                            segmentTimestamps(currSegment + 1) = timestamps(vals(jj) + 1);
                            currSegment = currSegment + 1;
                        end
                        
                        % restore file pointer position
                        fseek(FID, currBytePosition, 'bof');
                    end
                    
                    % update for next round
                    % -1 on the number of packets processed because the last packet
                    % is included in the next frame also
                    timestampFirst = timestampLast;
                    numPacketsProcessed = numPacketsProcessed + frameNumPackets - 1;
                end
                numPacketsProcessed = numPacketsProcessed + 1; % account for the overlapped sample on each frame
                assert(numPacketsProcessed == numPacketsTotal, 'Inconsistent number of packets processed (%d) versus number of packets in file (%d)',numPacketsProcessed,(f.EOF-f.EOexH)/packetSize);
                
                % compute number of datapoints in the last segment
                % add one to the number of packets processed to account for the
                % last packet of the file not being included in a subsequent frame
                segmentDatapoints(currSegment) = numPacketsProcessed - sum(segmentDatapoints(~isnan(segmentDatapoints)));
                segmentDurations(currSegment) = timestampLast - segmentTimestamps(currSegment) + 1;
    
                % add into NSx structure
                NSx.MetaTags.Timestamp = segmentTimestamps(1:currSegment);
                NSx.MetaTags.DataPoints = segmentDatapoints(1:currSegment);
                NSx.MetaTags.DataDurationSec = segmentDurations(1:currSegment)/NSx.MetaTags.TimeRes;
                file.MetaTags.DataDurationTimeRes = segmentDurations(1:currSegment);
            else

                % add into NSx structure
                fseek(FID,f.EOexH + 1,'bof'); % + byte (header)
                NSx.MetaTags.Timestamp = fread(FID,1,timestampType);
                NSx.MetaTags.DataPoints = (f.EOF - f.EOexH)/packetSize;
                NSx.MetaTags.DataDurationSec = NSx.MetaTags.DataPoints/NSx.MetaTags.SamplingFreq;
                file.MetaTags.DataDurationTimeRes = NSx.MetaTags.DataPoints*NSx.MetaTags.TimeRes/NSx.MetaTags.SamplingFreq;
            end
        else
            segmentCount = 0;
            while double(ftell(FID)) < f.EOF
                headerByte = fread(FID, 1, 'uint8=>double');
                if headerByte ~= 1
                    % Fixing another bug in Central 6.01.00.00 TOC where DataPoints is
                    % not written back into the Data Header
                    %% BIG NEEDS TO BE FIXED
                    NSx.MetaTags.DataPoints = floor(double(f.EOF - (f.EOexH+1+timestampSize+4))/(channelCount*2));
                    NSx.MetaTags.DataDurationSec = NSx.MetaTags.DataPoints/NSx.MetaTags.SamplingFreq;
                    break;
                end
                segmentCount = segmentCount + 1;
                startTimestamp = fread(FID, 1, timestampType);
                if flagMultiNSP
                    
                    % close existing (read-only) file descriptor
                    currBytePosition = ftell(FID);
                    fclose(FID);
                    
                    % open a file descriptor for read/write, write, and close
                    FID = fopen([requestedFilePath requestedFileName], 'r+', 'ieee-le');
                    startTimestamp = startTimestamp + syncShift;
                    fseek(FID, -timestampSize, 'cof');
                    fwrite(FID, startTimestamp, '*uint32');
                    fclose(FID);
                    
                    % re-open read-only and seek to remembered position
                    FID = fopen([requestedFilePath requestedFileName], 'r', 'ieee-le');
                    fseek(FID,currBytePosition,'bof');
                end
                NSx.MetaTags.Timestamp(segmentCount) = startTimestamp;
                NSx.MetaTags.DataPoints(segmentCount) = fread(FID, 1, 'uint32=>double');
                NSx.MetaTags.DataDurationSec(segmentCount) = NSx.MetaTags.DataPoints(segmentCount)/NSx.MetaTags.SamplingFreq;
                file.MetaTags.DataDurationTimeRes(segmentCount) = NSx.MetaTags.DataPoints(segmentCount)*NSx.MetaTags.TimeRes/NSx.MetaTags.SamplingFreq;
                fseek(FID, NSx.MetaTags.DataPoints(segmentCount) * channelCount * 2, 'cof');

                % Fixing the bug in 6.01.00.00 TOC where DataPoints is not
                % updated and is left as 0
                % NSx.MetaTags.DataPoints(segmentCount) = (f.EOData(segmentCount)-f.BOData(segmentCount))/(ChannelCount*2);
            end
        end
    end
    
    %% Calculate file pointers for data
    if strcmpi(NSx.MetaTags.FileTypeID, 'NEURALSG')
        % Determining DataPoints
        f.BOData = f.EOexH;
        f.EOData = f.EOF;
    elseif or(strcmpi(NSx.MetaTags.FileTypeID, 'NEURALCD'), strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP'))
        byteOffset = 1 + timestampSize + 4;
        if flagOneSamplePerPacket
            segmentOffset = f.EOexH;
            f.BOData = segmentOffset + byteOffset + [0 packetSize*cumsum(NSx.MetaTags.DataPoints(1:end-1))];
            f.EOData = segmentOffset + packetSize*cumsum(NSx.MetaTags.DataPoints);
        else
            segmentOffset = f.EOexH + (1:length(NSx.MetaTags.DataPoints))*byteOffset;
            f.BOData = segmentOffset + [0 cumsum(channelCount*NSx.MetaTags.DataPoints(1:end-1)*2)];
            f.EOData = segmentOffset + 2*channelCount*cumsum(NSx.MetaTags.DataPoints) - 1;
        end
    end
    
    % Determining if the file has a pause in it
    if length(NSx.MetaTags.DataPoints) > 1
        NSx.RawData.PausedFile = 1;
    end
    
    %% Save data headers for saveNSx
    dataHeaderBytes = cell(1,length(f.BOData));
    for ss = 1:length(f.BOData)
        headerByteSize = 1 + timestampSize + 4;
        fseek(FID, f.BOData(ss)-headerByteSize, 'bof');
        dataHeaderBytes{ss} = fread(FID, headerByteSize, '*uint8');
    end
    NSx.RawData.DataHeader = cat(1,dataHeaderBytes{:});
    
    %% Identify and validate requested first and last data points
    % Note that new files that use precision time protocol (PTP), for which
    % there is one sample per packet and hence the user request translates
    % to a starting and ending data packet. In older NSx files, there were
    % multiple samples per packet and the user request is for a starting
    % sample (which would exist inside a data packet). This code refers to
    % "data points" as a broad term to encapsulate both scenarios.
    if ~flagModifiedTime
        
        % default whole file
        requestedStartDataPoint = 1;
        requestedEndDataPoint = sum(NSx.MetaTags.DataPoints);
    else
        
        % TO-DO: utilize Gemini sample-by-sample timestamps
        switch requestedTimeScale
            case {'sec', 'secs', 'second', 'seconds'}
                
                % convert seconds to samples
                requestedStartDataPoint = requestedStartValue * NSx.MetaTags.SamplingFreq + 1;
                requestedEndDataPoint = requestedEndValue * NSx.MetaTags.SamplingFreq;
            case {'min', 'mins', 'minute', 'minutes'}
                
                % convert minutes to samples
                requestedStartDataPoint = requestedStartValue * NSx.MetaTags.SamplingFreq * 60 + 1;
                requestedEndDataPoint = requestedEndValue * NSx.MetaTags.SamplingFreq * 60;
            case {'hour', 'hours'}
                
                % convert hours to samples
                requestedStartDataPoint = requestedStartValue * NSx.MetaTags.SamplingFreq * 3600 + 1;
                requestedEndDataPoint = requestedEndValue * NSx.MetaTags.SamplingFreq * 3600;
            case {'sample','samples'}
                
                % carry over samples
                requestedStartDataPoint = requestedStartValue;
                requestedEndDataPoint = requestedEndValue;
            otherwise
                
                % should never get here based on input argument processing
                error('Unknown requested time scale');
        end
    end
    requestedNumDataPoints = requestedEndDataPoint - requestedStartDataPoint + 1;
    
    % validate start and end data points
    assert(requestedEndDataPoint>=requestedStartDataPoint,'Start data point (%d) must be less than the end data point (%d)',requestedStartDataPoint,requestedEndDataPoint);
    assert(requestedStartDataPoint<=sum(NSx.MetaTags.DataPoints),'Start data point (%d) is greater than total number of data point (%d)',requestedStartDataPoint,sum(NSx.MetaTags.DataPoints));
    if requestedStartDataPoint <= 0
        warning('Start data point (%d) must be greater than or equal to 1; updating to comply.',requestedStartDataPoint);
        requestedStartDataPoint = 1;
    end
    if requestedEndDataPoint > sum(NSx.MetaTags.DataPoints)
        warning('End data point (%d) must be less than or equal to the total number of data points (%d).',requestedEndDataPoint,sum(NSx.MetaTags.DataPoints));
        response = input('Do you wish to update the last requested data point to last one available in the file and continue? (y/N) ', 's');
        if strcmpi(response,'y')
            warning('Changed end data point from %d to %d',requestedEndDataPoint,sum(NSx.MetaTags.DataPoints));
            requestedEndDataPoint = sum(NSx.MetaTags.DataPoints);
        else
            error('Invalid last requested data point');
        end
    end
    
    %% Identify data segments containing requested packets
    requestedSegments = nan(1,2); % first and last requested segments
    startTimeStampShift = 0;
    
    % user requested specific data points: look for start/end segments
    % containing these data points
    if flagModifiedTime

        % search first for segment containing the requested data points
        dataPointOfInterest = requestedStartDataPoint;
        segmentStartDataPoint = zeros(1,length(NSx.MetaTags.DataPoints));
        segmentDataPoints = zeros(1,length(NSx.MetaTags.DataPoints));
        
        % loop over data segments
        for currSegment = 1:length(NSx.MetaTags.DataPoints)
            if dataPointOfInterest <= sum(NSx.MetaTags.DataPoints(1:currSegment))

                % still looking at the starting value
                if all(isnan(requestedSegments))

                    % set the starting data point
                    segmentStartDataPoint(currSegment) = dataPointOfInterest;

                    % check if end point also in this segment
                    if requestedEndDataPoint <= sum(NSx.MetaTags.DataPoints(1:currSegment))

                        % set number of data points read for this segment
                        segmentDataPoints(currSegment) = requestedEndDataPoint - segmentStartDataPoint(currSegment) + 1;

                        % set the start/end segments to this one
                        requestedSegments = [currSegment currSegment];

                        % exit out of loop
                        break;
                    end

                    % since end point not here, read all remaining points
                    segmentDataPoints(currSegment) = NSx.MetaTags.DataPoints(currSegment) - dataPointOfInterest + 1;

                    % switch to looking for the end point
                    dataPointOfInterest = requestedEndDataPoint;

                    % compute shift applied to timestamp for this segment
                    startTimeStampShift = (segmentStartDataPoint(currSegment)-1) * NSx.MetaTags.TimeRes / NSx.MetaTags.SamplingFreq - NSx.MetaTags.Timestamp(currSegment);
                else

                    % segment containing the ending data point
                    segmentStartDataPoint(currSegment) = 1;
                    segmentDataPoints(currSegment) = requestedNumDataPoints - sum(segmentDataPoints(1:currSegment-1));
                    requestedSegments(2) = currSegment;
                    break;
                end

                % save out current segment as containing the start point
                requestedSegments(1) = currSegment;
            else

                % this segment does not contain the point of interest
                if all(isnan(requestedSegments))

                    % haven't found start or end segment yet
                    segmentStartDataPoint(currSegment) = NSx.MetaTags.DataPoints(currSegment);
                    segmentDataPoints(currSegment) = 0;
                elseif isnan(requestedSegments(2))

                    % already found the start segment, but not end
                    segmentStartDataPoint(currSegment) = 1;
                    segmentDataPoints(currSegment) = NSx.MetaTags.DataPoints(currSegment);
                else

                    % found both start and end already
                    segmentStartDataPoint(currSegment) = 1;
                    segmentDataPoints(currSegment) = 0;
                end
            end
        end
    else

        % no user input; read the whole file
        requestedSegments = [1 length(NSx.MetaTags.DataPoints)];
        segmentStartDataPoint = ones(1,length(NSx.MetaTags.DataPoints));
        segmentDataPoints = NSx.MetaTags.DataPoints;
    end
    
    % calculate total number of data points requested
    numDataPointsRequested = sum(floor(segmentDataPoints/requestedSkipFactor));
    
    %% Read data
    file.MetaTags.DataPoints = NSx.MetaTags.DataPoints;
    file.MetaTags.DataDurationSec = NSx.MetaTags.DataDurationSec;
    file.MetaTags.Timestamp = NSx.MetaTags.Timestamp;
    numDataPointsRead = 0;
    if flagReadData
        
        % loop over requested data segments
        NSx.Data = cell(1,diff(requestedSegments)+1);
        for currSegment = requestedSegments(1):requestedSegments(2)
            outputSegment = currSegment - requestedSegments(1) + 1;
            
            % seek to start of data
            fseek(FID, f.BOData(currSegment), 'bof');
            
            % seek to first requested packet in the current segment
            if flagOneSamplePerPacket
                fseek(FID, (segmentStartDataPoint(currSegment) - 1) * packetSize, 'cof');
            else
                fseek(FID, (segmentStartDataPoint(currSegment) - 1) * 2 * channelCount, 'cof');
            end
            
            % seek to first requested channel in the current packet
            fseek(FID, (requestedFirstChannel-1) * 2, 'cof');
            
            % set up parameters for reading data
            precisionString = sprintf('%d*int16=>%s',numChannelsToRead,requestedPrecisionType);
            outputDimensions = [numChannelsToRead floor(segmentDataPoints(currSegment)/requestedSkipFactor)];
            if flagOneSamplePerPacket
                bytesToSkipNormal = packetSize - 2*numChannelsToRead; % standard (i.e., skip factor==1)
                bytesSkipFactor = packetSize*(requestedSkipFactor - 1); % additional to skip (skip factor > 1)
            else
                bytesToSkipNormal = 2*(channelCount - numChannelsToRead);
                bytesSkipFactor = 2*channelCount*(requestedSkipFactor-1);
            end
            bytesToSkip = bytesToSkipNormal + bytesSkipFactor; % total
            
            % read data
            NSx.Data{outputSegment} = fread(FID, outputDimensions, precisionString, bytesToSkip);
        end

        % read timestamps - loop over the requested data segments
        if (~flagAlign || ~flagSegment) && flagOneSamplePerPacket
            NSx.Time = cell(1,diff(requestedSegments)+1);
            for currSegment = requestedSegments(1):requestedSegments(2)
                outputSegment = currSegment - requestedSegments(1) + 1;

                % seek to start of data
                fseek(FID, f.BOData(currSegment), 'bof');

                % seek to first requested packet in the current segment
                fseek(FID, (segmentStartDataPoint(currSegment) - 1) * packetSize, 'cof');

                % seek back to the timestamp
                fseek(FID, -(4 + timestampSize), 'cof');

                % set up parameters for reading data
                precisionString = sprintf('%s=>double',timestampType);
                outputDimensions = [1 floor(segmentDataPoints(currSegment)/requestedSkipFactor)];
                bytesToSkipNormal = packetSize - timestampSize; % standard (i.e., skip factor==1)
                bytesSkipFactor = packetSize*(requestedSkipFactor - 1); % additional to skip (skip factor > 1)
                bytesToSkip = bytesToSkipNormal + bytesSkipFactor; % total

                % read data
                NSx.Time{outputSegment} = fread(FID, outputDimensions, precisionString, bytesToSkip);
            end
        end

        % define user tags: info specific to data being read
        NSx.MetaTags.Timestamp = NSx.MetaTags.Timestamp(requestedSegments(1):requestedSegments(2));
        NSx.MetaTags.Timestamp(1) = NSx.MetaTags.Timestamp(1) + startTimeStampShift;
        NSx.MetaTags.DataPoints = cellfun(@(x) size(x,2), NSx.Data, 'UniformOutput', true);
        NSx.MetaTags.DataDurationSec = NSx.MetaTags.DataDurationSec(requestedSegments(1):requestedSegments(2));
        NSx.MetaTags.DataDurationSec(1) = NSx.MetaTags.DataDurationSec(1) - (segmentStartDataPoint(requestedSegments(1))-1)/NSx.MetaTags.SamplingFreq;
        NSx.MetaTags.DataDurationSec(end) = NSx.MetaTags.DataPoints(end)/NSx.MetaTags.SamplingFreq;

        % verify amount of data read from disk
        numDataPointsRead = sum(cellfun(@(x)size(x,2),NSx.Data));
        if numDataPointsRead~=numDataPointsRequested
            warning('Expected to read %d data points, but output has %d data points.',numDataPointsRequested,numDataPointsRead);
        end
    end
catch ME
    fclose(FID);
    rethrow(ME);
end
fclose(FID);

%% Bug fix
% Fix a bug in 6.03 where data packets with 0 length may be added
if flagReadData && any(NSx.MetaTags.DataPoints == 0)
    segmentsThatAreZero = find(NSx.MetaTags.DataPoints == 0);
    NSx.MetaTags.DataPoints(segmentsThatAreZero) = [];
    NSx.MetaTags.DataDurationSec(segmentsThatAreZero) = [];
    NSx.MetaTags.Timestamp(segmentsThatAreZero) = [];
    NSx.Data(segmentsThatAreZero) = [];
    if isfield(NSx,'Time')
        NSx.Time(segmentsThatAreZero) = [];
    end
end

%% Remove extra channels that were read, but weren't supposed to be read
if flagReadData
    channelsRead = min(requestedChannelIndex):max(requestedChannelIndex);
    idxToKeep = ismember(channelsRead,requestedChannelIndex);
    NSx.Data = cellfun(@(x)x(idxToKeep,:),NSx.Data,'UniformOutput',false);
end
if isfield(NSx,'ElectrodesInfo')
    idxToRemove = ~ismember(1:length(NSx.ElectrodesInfo),requestedChannelIndex);
    NSx.ElectrodesInfo(idxToRemove) = [];
end
if isfield(NSx.MetaTags,'ChannelID')
    idxToRemove = ~ismember(1:length(NSx.MetaTags.ChannelID),requestedChannelIndex);
    NSx.MetaTags.ChannelID(idxToRemove) = [];
end

%% Zero-pad data if requested
if flagReadData && flagZeroPad
    
    % only operate on first data segment
    currSegment = 1;
    
    % compute how many zeros and total number of values to add
    numZerosToAdd = floor(NSx.MetaTags.Timestamp(currSegment) / requestedSkipFactor);
    if flagMultiNSP
        numZerosToAdd = numZerosToAdd + syncShift;
    end
    numValuesToAdd = NSx.MetaTags.ChannelCount * numZerosToAdd;
    
    % sanity check
    if strcmpi(NSx.MetaTags.FileTypeID, 'BRSMPGRP')
        
        % calculate how many values and how many bytes
        numValuesOfData = sum(cellfun(@numel,NSx.Data));
        if strcmpi(requestedPrecisionType,'int16')
            numBytesToAdd = numValuesToAdd*2;
            numBytesOfData = numValuesOfData*2;
        elseif strcmpi(requestedPrecisionType,'double')
            numBytesToAdd = numValuesToAdd*8;
            numBytesOfData = numValuesOfData*8;
        end
        
        % check whether to show the warning
        flagShowZeroPadWarning = 1;
        if flagFoundSettingsManager
            flagShowZeroPadWarning = NPMKSettings.ShowZeroPadWarning;
        end
        
        % generate warning and ask to continue
        if NSx.MetaTags.Timestamp(1) > 30000 && flagShowZeroPadWarning
            warning('Zero padding would add %d bytes to data which already require %d bytes in memory', numBytesToAdd, numBytesOfData);
            response = input('Do you wish to continue? (y/N) ', 's');
            if ~strcmpi(response, 'y')
                warning('Turned off zero padding by user request');
                flagZeroPad = 0;
            end
        end
        
        % check on continuing to show this warning
        if flagFoundSettingsManager
            response = input('Do you want NPMK to continue to ask you about this every time? (Y/n) ', 's');
            if ~strcmpi(response, 'n')
                NPMKSettings.ShowZeroPadWarning = 1;
            else
                NPMKSettings.ShowZeroPadWarning = 0;
            end
            settingsManager(NPMKSettings);
        end
    end
    
    % perform zero padding
    % (flagZeroPad may be set to false above, so need to re-evaluate)
    if requestedStartDataPoint == 1 && flagZeroPad
        
        % only for the first data segment
        NSx.Data{currSegment} = [zeros(NSx.MetaTags.ChannelCount, numZerosToAdd, requestedPrecisionType) NSx.Data{currSegment}];
        if isfield(NSx,'Time')
            NSx.Time{currSegment} = [zeros(1, numZerosToAdd, 'double') NSx.Time{currSegment}];
        end

        % update metadata
        NSx.MetaTags.DataDurationSec(currSegment) = size(NSx.Data{currSegment},2)/NSx.MetaTags.SamplingFreq;
        NSx.MetaTags.DataPoints(currSegment) = size(NSx.Data{currSegment},2);
    end
    NSx.MetaTags.Timestamp(currSegment) = 0;
end


%% Adjust for the data's unit.
if flagReadData 
    if flagConvertToUv
        NSx.Data = cellfun(@(x) bsxfun(@rdivide, x, 1./(double([NSx.ElectrodesInfo.MaxAnalogValue])./double([NSx.ElectrodesInfo.MaxDigiValue]))'),NSx.Data ,'UniformOutput',false);
    else
        flagShowuVWarning = 1;
        if flagFoundSettingsManager
            flagShowuVWarning = NPMKSettings.ShowuVWarning;
        end
        if flagShowuVWarning
            warning('Note that data are in units of 1/4 µV; see ''uv'' argument.');
        end
        if flagShowuVWarning && flagFoundSettingsManager
            response = input('Do you want NPMK to continue to warn you about this every time? (Y/n) ', 's');
            if ~strcmpi(response, 'n')
                NPMKSettings.ShowuVWarning = 1;
            else
                NPMKSettings.ShowuVWarning = 0;
            end
            settingsManager(NPMKSettings);
        end
    end
end

%% Add implementation of samplealign for cases where it is needed
if flagReadData && flagOneSamplePerPacket
    if flagAlign
        for ii = 1:length(NSx.Data)
            fileDataLength = file.MetaTags.DataPoints(ii);
            fileDuration = file.MetaTags.DataDurationTimeRes(ii);
    
            % Calculate the ratio between time gaps and expected time gap
            % based on the sampling rate of the recording. A recording
            % where the claimed sampling rate and true sampling rate based
            % off PTP time are identical will have a ratio of 1;
            samplingRates = fileDuration/fileDataLength/NSx.MetaTags.TimeRes*NSx.MetaTags.SamplingFreq;
    
            % Calculate the number of samples that should be added or
            % removed
            addedSamples = round((samplingRates-1)*fileDataLength);
    
            % Establish where the points should be added or removed
            gapIndex = round(fileDataLength/(abs(addedSamples)+1));
    
            % calculate the portion of samples added/subtracted to the
            % requested data, which may be shorter than the full file
            % use floor because we need addedsamples+1 sections to avoid
            % adding/subtracting samples at the beginning or end of the data.
            addedSamples = floor(addedSamples * NSx.MetaTags.DataPoints(ii)/fileDataLength);
            if addedSamples == 0
                continue;
            end
            
            % split into cell arrays
            dim1Size = size(NSx.Data{ii},1);
            if gapIndex >= size(NSx.Data{ii},2)
                if abs(addedSamples)>1
                    warning('Expected to add or remove only one sample')
                end
                dim2Size = [round(size(NSx.Data{ii},2)/2) size(NSx.Data{ii},2)-round(size(NSx.Data{ii},2)/2)];
            else
                dim2Size = [repmat(gapIndex,1,abs(addedSamples)) size(NSx.Data{ii},2) - gapIndex*abs(addedSamples)];
            end
            NSx.Data{ii} = mat2cell(NSx.Data{ii},dim1Size,dim2Size);
    
            % add or subtract
            if abs(addedSamples)==1
                sampleString = sprintf('%d sample',abs(addedSamples));
                whereString = 'at midpoint';
            else
                sampleString = sprintf('%d samples',abs(addedSamples));
                whereString = 'evenly spaced';
            end
            if length(NSx.Data)==1
                segmentString = 'the data';
            else
                segmentString = sprintf('data segment %d/%d',ii,length(NSx.Data));
            end
            if addedSamples>0
                NSx.Data{ii}(1:end-1) = cellfun(@(x) [x x(:,end)], NSx.Data{ii}(1:end-1), 'UniformOutput',false);
                warning('Added %s to %s (%s) for clock drift alignment',sampleString,segmentString,whereString)
            elseif addedSamples<0
                NSx.Data{ii}(1:end-1) = cellfun(@(x) x(:,1:end-1), NSx.Data{ii}(1:end-1), 'UniformOutput',false);
                warning('Removed %s from %s (%s) for clock drift alignment',sampleString,segmentString,whereString)
            end
    
            % combine to form the full data again
            NSx.Data{ii} = cat(2,NSx.Data{ii}{:});
    
            % recompute some metadata
            NSx.MetaTags.DataPoints(ii) = size(NSx.Data{ii},2);
            NSx.MetaTags.DataDurationSec(ii) = size(NSx.Data{ii},2)/NSx.MetaTags.SamplingFreq;
        end
    end
end

% remove Time field when not used
if ~flagOneSamplePerPacket || (flagSegment && flagAlign)
    NSx = rmfield(NSx,'Time');
end

% reduce to array if only one cell
if flagReadData && iscell(NSx.Data) && length(NSx.Data)==1
    NSx.Data = NSx.Data{1};
    if isfield(NSx,'Time')
        NSx.Time = NSx.Time{1};
    end
end

% Display a report of basic file information and the Basic Header.
if flagReport
    disp( '*** FILE INFO **************************');
    disp(['File Path          = '  NSx.MetaTags.FilePath]);
    disp(['File Name          = '  NSx.MetaTags.Filename]);
    disp(['File Extension     = '  NSx.MetaTags.FileExt]);
    disp(['File Version       = '  NSx.MetaTags.FileSpec]);
    disp(['Duration (seconds) = '  num2str(NSx.MetaTags.DataDurationSec)]);
    disp(['Total Datapoints   = '  num2str(NSx.MetaTags.DataPoints)]);
    disp(' ');
    disp( '*** BASIC HEADER ***********************');
    disp(['File Type ID       = '  NSx.MetaTags.FileTypeID]);
    disp(['Sample Frequency   = '  num2str(double(NSx.MetaTags.SamplingFreq))]);
    disp(['Electrodes Read    = '  num2str(double(NSx.MetaTags.ChannelCount))]);
    disp(['Datapoints Read    = '  num2str(numDataPointsRead)]);
end

% Create output variable in user workspace even if no output argument
outputName = ['NS' requestedFileExtension(4)];
if (nargout == 0)
    assignin('caller', outputName, NSx);
else
    varargout{1} = NSx;
end

% Print the load time
if flagReport
    disp(['The load time for ' outputName ' file was ' num2str(toc, '%0.1f') ' seconds.']);
end

end