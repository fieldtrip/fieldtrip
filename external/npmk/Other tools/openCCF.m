function [infoPackets, version] = openCCF(filename)

% openCCF
%
% This script reads a .ccf file and outputs it in a structure.
%
%   filename:  Name of the file to be opened. If the fname is omitted
%              the user will be prompted to select a file. 
%              DEFAULT: Will open Open File UI.
%
%   Channels   1-128 are neural channels
%   Channels 129-144 are external analog inputs
%   Channels 145-148 are external analog outpus
%   Channels
%    
%   Kian Torab
%   kian@blackrockmicro.com
%   Blackrock Microsystems
%   Version 2.1.0.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
% 1.2.3.0:
%   - Minor bug fix that led to a crash in certain cases.
%
% 1.2.4.0:
%   - Minor bug fix regarding passing a filename variable to the function.
%
% 2.0.0.0: 
%   - Implemented XML CCF file format.
%
% 2.1.0.0:
%   - Fixed a bug in loading nTrode groups with a base of 0.
%
% 2.2.0.0 April 29, 2020
%   - Fixed an error where N-Trodes with less than 4 members read an extra
%     1 as the extra non-existent members.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infoPackets.openCCFVersion = '2.1.0.0';

%% Parameters

% Which versions of .ccf's are supported
supportedVersions = {'3.6', '3.7', '3.8', '3.9'};
% Number of cbPKT_CHANINFO packets to read in the file
nValidPKT_CHANINFO = 160;

if( nargin ~= 1 )
    [filename, pathname] = getFile('*.ccf');
    if ~filename
        disp('No file was selected.');
        return;
    end
    fullfilename = [pathname filename];
elseif (nargin == 1)
    fullfilename = filename;
else
    if strcmpi(filename, 'ver')
    disp(['Version is ' infoPackets.openCCFVersion]);
    infoPackets = [];
    return;
    end
end


%% Open file
[fid, message] = fopen( fullfilename, 'rb', 'l');
if( fid == -1 ),
   warning(['Unable to open file: ' fullfilename '. Error message ' message]);
   return;
end;

%% Read header incl. version
head = fread(fid, 16, 'uint8=>char*1' )';

if strcmp(head(1:5), 'cbCCF')
    version = deblank(strtrim(head(6:end))); % Strip possible leading space and trailing nulls
elseif strcmp(head(1:5), '<?xml')
    version = '3.9';
else
  warning('Not a well-formed .ccf file');
  fclose(fid);
  return;
end

% Save version info in the structure.
infoPackets.Version = version;

% Verify to see if the file version is supported
if ~ismember(version, supportedVersions)
  warning(sprintf('Unsupported file version: %s', version));
  fclose(fid);
  return;
end

if strcmpi(version, '3.9')
    infoPackets = parseCCF(fullfilename);
    for nTrodeIDX = 1:length(infoPackets.Children(7).Children)
        for trodeIDX = 1:str2double(infoPackets.Children(7).Children(nTrodeIDX).Children(7).Children.Data)
            tempTrodes(trodeIDX) =  ...
            str2double(infoPackets.Children(7).Children(nTrodeIDX).Children(9).Children(trodeIDX).Children.Data) + 1;
        end
        if ~isempty(tempTrodes)
            infoPackets.NTrodeInfo.NTrodeID(nTrodeIDX) = nTrodeIDX;
            infoPackets.NTrodeInfo.NTrodeMembers{nTrodeIDX} = tempTrodes;
        end 
        tempTrodes = [];
    end
    return;
end

%% Read packets
for p = 1:nValidPKT_CHANINFO
    if p == 128
%    disp('w');
end

  infoPackets.ChanInfo(p) = readInfoPacket(fid, version);
end
infoPackets.AdaptiveFilterInfo = readAdapFilterPacket(fid);
infoPackets.SortingInfo = readSortingPacket(fid);
infoPackets.SystemInfo = []; %readSystemInfoPacket(fid);
infoPackets.NTrodeInfo = readNTrodeInfoPacket(infoPackets);

%% Close file
fclose(fid);



function p = readInfoPacket(fid, version)
% Read a cbPKT_CHANINFO packet
    p = [];
    cbLEN_STR_LABEL = 16;
    cbMAXUNITS = 5;
    cbMAXHOOPS = 4;
    p.SystemClockTimestamp = fread(fid, 1, '*uint32');
    p.Child = fread(fid, 1, '*uint16');
    p.AnaInPacketType = fread(fid, 1, '*uint8');
    p.SystemDataLength = fread(fid, 1, '*uint8');
    p.AnaInChannelID = fread(fid, 1, '*uint32');
    p.AnaInChannelProcessor = fread(fid, 1, '*uint32');
    p.AnaInBankID = fread(fid, 1, '*uint32');
    p.AnaInPinID = fread(fid, 1, '*uint32');
    h = fread(fid, 1, '*uint32');
    p.CapabilitiesGeneralChannel = setFlagsCBCHAN(h);
    h = fread(fid, 1, '*uint32');
    p.CapabilitiesDigOutChannel = setFlagsCBDOUT(h);
    h = fread(fid, 1, '*uint32');
    p.CapabilitiesDigInChannel = setFlagsCBDINP(h);
    h = fread(fid, 1, '*uint32');
    p.CapabilitiesAnaOutChannel = setFlagsCBAOUT(h);
    h = fread(fid, 1, '*uint32');
    p.CapabilitiesAnaInChannel = setFlagsCBAINP(h);
    p.CapabilitiesSpikeProcessing = fread(fid, 1, '*uint32');

    p.InputPhysicalChannelScalingInfo       = readCbScaling(fid);
    p.InputPhysicalChannelFilterDefinition  = readFiltDesc(fid);
    p.OutputPhysicalChannelScalingInfo      = readCbScaling(fid);
    p.OutputPhysicalChannelFilterDefinition = readFiltDesc(fid);

    p.ChannelLabel     = fread(fid, cbLEN_STR_LABEL, '*char*1')';
    p.ChannelUserFlags = fread(fid, 1, '*uint32');
    p.ChannelPosition  = fread(fid, 4, '*int32')';

    p.InputPhysicalChannelScalingInfoCustom  = readCbScaling(fid);
    p.OuputPhysicalChannelScalingInfoCustom = readCbScaling(fid);

    p.doutopts = fread(fid, 1, '*uint32');
    p.dinpopts = fread(fid, 1, '*uint32');
    p.aoutopts = fread(fid, 1, '*uint32');
    p.eopchar  = fread(fid, 1, '*uint32');

    % Below two lines are actually part of a union, not sure how to figure out
    % which option is actually being used
    p.AnaInMonitorChannelAddress             = fread(fid, 1, '*uint32');
    p.AnaInOutputValue                       = fread(fid, 1, '*int32');
    p.AnaInOptions                           = fread(fid, 1, '*uint32');
    p.AnaInLNCAdaptationRate                 = fread(fid, 1, '*uint32');
    p.AnaInContinuousStreamFilterID          = fread(fid, 1, '*uint32');
    p.AnalogInputContinuousStreamSampleGroup = fread(fid, 1, '*uint32');
    p.AnaInContinuousStreamDisplayFactorMin  = fread(fid, 1, '*int32');
    p.AnaInContinuousStreamDisplayFactorMax  = fread(fid, 1, '*int32');
    p.AnaInSpikeStreamFilterID               = fread(fid, 1, '*uint32');
    p.AnaInSpikeStreamDisplayFactorMax       = fread(fid, 1, '*int32');
    p.AnaInLNCDisplayFactorMax               = fread(fid, 1, '*int32');
    p.AnaInSpikeProcessingOptions            = fread(fid, 1, '*uint32');
    p.AnaInSpikeStreamThresholdLevel         = fread(fid, 1, '*int32');
    p.AnaInSpikeStreamThresholdLimit         = fread(fid, 1, '*int32');
    p.AnaInNTrodeGroupID                     = fread(fid, 1, '*uint32');
    p.AnaInAmplitudeRejectionValuePositive   = fread(fid, 1, '*int16');
    p.AnaInAmplitudeRejectionValueNegative   = fread(fid, 1, '*int16');
    p.AnaInDigitalReferencingChannel         = fread(fid, 1, '*uint32');
    p.AnaInManualMappingUnit                 = readCbManualUnitMappings(fid, cbMAXUNITS, version);
    p.AnaInSpikeHoopsSorting                 = readCbHoops(fid, cbMAXUNITS, cbMAXHOOPS);

function p = readAdapFilterPacket(fid)
    p = [];
    p.Time = fread(fid, 1, '*uint32');
    p.Child = fread(fid, 1, '*uint16');
    p.Type = fread(fid, 1, '*uint8');
    p.DataLength = fread(fid, 1, '*uint8');
    p.Chan = fread(fid, 1, '*uint32');
    p.Mode = fread(fid, 1, '*uint32');
    p.LearningRate = fread(fid, 1, '*float32');
    p.RefChan1 = fread(fid, 1, '*uint32');
    p.RefChan2 = fread(fid, 1, '*uint32');


function p = readSortingPacket(fid)
    p = [];
    p.Time = fread(fid, 1, '*uint32');
    p.Child = fread(fid, 1, '*uint16');
    p.Type = fread(fid, 1, '*uint8');
    p.DataLength = fread(fid, 1, '*uint8');
    p.NumMaxSimultChans = fread(fid, 1, '*uint32');
    p.RefractoryPeriodInSamples = fread(fid, 1, '*uint32');


function p = readSystemInfoPacket(fid)
    p = [];

function p = readNTrodeInfoPacket(infoPackets)
    p = [];
    NTrodeGroups = unique([infoPackets.ChanInfo.AnaInNTrodeGroupID]); NTrodeGroups(1) = [];
    NTrodeGroupsCount = length(NTrodeGroups);

for ntrodeIDX = 1:NTrodeGroupsCount
    p.NTrodeID(ntrodeIDX) = NTrodeGroups(ntrodeIDX);
    p.NTrodeMembers{ntrodeIDX} = find([infoPackets.ChanInfo.AnaInNTrodeGroupID] == ntrodeIDX);
end

%% Sub-struct readers
function sc = readCbScaling(fid)
% Read a cbSCALING struct

cbLEN_STR_UNIT = 8;

sc.AnalogChannelScaleMin = fread(fid, 1, '*int16');
sc.AnalogChannelScaleMax = fread(fid, 1, '*int16');
sc.AnalogChannelScaleMin = fread(fid, 1, '*int32');
sc.AnalogChannelScaleMax = fread(fid, 1, '*int32');
sc.AnalogChannelGain = fread(fid, 1, '*int32');
sc.AnalogChannelUnit = fread(fid, cbLEN_STR_UNIT, '*char*1')';



function fd = readFiltDesc(fid)
% Read a cbFILTDESC struct

cbLEN_STR_FILT_LABEL = 16;

fd.FilterLabel = fread(fid, cbLEN_STR_FILT_LABEL, '*char*1')';
fd.FilterHighPassFrequency = fread(fid, 1, '*uint32');
fd.FilterHighPassOrder = fread(fid, 1, '*uint32');
fd.FilterHighPassType = fread(fid, 1, '*uint32');
fd.FilterLowPassFrequency = fread(fid, 1, '*uint32');
fd.FilterLowPassOrder = fread(fid, 1, '*uint32');
fd.FilterLowPassType = fread(fid, 1, '*uint32');



function ums = readCbManualUnitMappings(fid, maxUnits, version)
% Loop to read all the cbMANUALUNITMAPPING packets. Note that this requires
% knowing the file version: v3.6 is different from 3.7 and 3.8

% Loop growing is inefficient, but it's just not that much data.
if strcmp(version, '3.6')
    for i = 1:maxUnits
        ums(i) = readCbManualUnitMapping3_6(fid);
    end
else
    for i = 1:maxUnits
        ums(i) = readCbManualUnitMapping(fid);
    end
end  


function um = readCbManualUnitMapping(fid)
% Read a cbMANUALUNITMAPPING packet for a v3.7 or 3.8 file

um.nOverride = fread(fid, 1,     '*int16');
um.afOrigin  = fread(fid, 3,     '*int16')';
um.afShape   = fread(fid, [3 3], '*int16')';
um.aPhi      = fread(fid, 1,     '*int16')';
um.bValid    = fread(fid, 1,     '*uint32')';

function um = readCbManualUnitMapping3_6(fid)
% Read a cbMANUALUNITMAPPING packet for a v3.6 file

um.nOverride = fread(fid, 1,     '*uint32');
um.afOrigin  = fread(fid, 3,     '*float32')';
um.afShape   = fread(fid, [3 3], '*float32')';
um.aPhi      = fread(fid, 1,     '*float32')';
um.bValid    = fread(fid, 1,     '*uint32')';




function hs = readCbHoops(fid, maxUnits, maxHoops)
% Loop to read the cbHOOP packets

% Loop growing is inefficient, but it's just not that much data.
for u = 1:maxUnits
  for h = 1:maxHoops
    hs(u, h) = readCbHoop(fid);
  end
end


function h = readCbHoop(fid)
% Read a cbHOOP packet
h.valid = fread(fid, 1, '*uint16')';
h.time = fread(fid, 1, '*int16')';
h.min = fread(fid, 1, '*int16')';
h.max = fread(fid, 1, '*int16')';

function h = setFlagsCBCHAN(flags)

binaryValue = dec2bin(flags,12);

% define constants
h.ChanExists     = 0;  % Channel id is allocated
h.ChanConnected  = 0;  % Channel is connected and mapped and ready to use
h.ChanIsolated   = 0;  % Channel is electrically isolated
h.ChanType       = '';

if strcmpi(binaryValue(1), '1')
    h.ChanType = 'DigitalOut';
end
if strcmpi(binaryValue(2), '1')
    h.ChanType = 'DigitalInput';
end
if strcmpi(binaryValue(3), '1')
    h.ChanType = 'AnalogOutput';
end
if strcmpi(binaryValue(4), '1')
    h.ChanType = 'AnalogInput';
end
if strcmpi(binaryValue(10), '1')
    h.ChanIsolated = 1;
end
if strcmpi(binaryValue(11), '1')
    h.ChanConnected = 1;
end
if strcmpi(binaryValue(12), '1')
    h.ChanExists = 1;
end

function h = setFlagsCBDOUT(flags)

binaryValue = dec2bin(flags,30);

h.BaudRate = 0;
h.NumberOfBits = 0;
h.CanBeManuallyConfigured = 0;
h.TrackMostRecentlySelected = 0;
h.CanOutputFrequency = 0;
h.ChanMonitorUnit0 = 0;
h.ChanMonitorUnit1 = 0;
h.ChanMonitorUnit2 = 0;
h.ChanMonitorUnit3 = 0;
h.ChanMonitorUnit4 = 0;
h.ChanMonitorUnit5 = 0;
h.ChanMonitorUnitAll = 0;

if strcmpi(binaryValue(1:6), '111111')
    h.ChanMonitorUnitAll = 1;
else
    if strcmpi(binaryValue(1), '1')
        h.ChanMonitorUnit5 = 1;
    end
    if strcmpi(binaryValue(2), '1')
        h.ChanMonitorUnit4 = 1;
    end
    if strcmpi(binaryValue(3), '1')
        h.ChanMonitorUnit3 = 1;
    end
    if strcmpi(binaryValue(4), '1')
        h.ChanMonitorUnit2 = 1;
    end
    if strcmpi(binaryValue(5), '1')
        h.ChanMonitorUnit1 = 1;
    end
    if strcmpi(binaryValue(6), '1')
        h.ChanMonitorUnit0 = 1;
    end
end
if strcmpi(binaryValue(12), '1')
    h.CanOutputFrequency = 1;
end
if strcmpi(binaryValue(13), '1')
    h.TrackMostRecentlySelected = 1;
end
if strcmpi(binaryValue(18), '1')
    h.CanBeManuallyConfigured = 1;
end
if strcmpi(binaryValue(19), '1')
    h.NumberOfBits = 32;
end
if strcmpi(binaryValue(20), '1')
    h.NumberOfBits = 16;
end
if strcmpi(binaryValue(21), '1')
    h.NumberOfBits = 8;
end
if strcmpi(binaryValue(22), '1')
    h.NumberOfBits = 1;
end
if strcmpi(binaryValue(25), '1')
    h.BaudRate = 115200;
end
if strcmpi(binaryValue(26), '1')
    h.BaudRate = 57600;
end
if strcmpi(binaryValue(27), '1')
    h.BaudRate = 38400;
end
if strcmpi(binaryValue(28), '1')
    h.BaudRate = 19200;
end
if strcmpi(binaryValue(29), '1')
    h.BaudRate = 9600;
end
if strcmpi(binaryValue(30), '1')
    h.BaudRate = 2400;
end

function h = setFlagsCBDINP(flags)

binaryValue = dec2bin(flags,22);

h.BaudRate = 0;
h.NumberOfBits = 0;
h.PortReadType = 0;
h.CapturePacketMethod = 0;
h.PortControlsEvent = 0;

if strcmpi(binaryValue(1), '1')
    h.PortReadType = '8-bit strobe/8-bit falling edge';
end
if strcmpi(binaryValue(2), '1')
    h.PortReadType = '8-bit strobe/8-bit rising edge';
end
if strcmpi(binaryValue(3), '1')
    h.PortReadType = '8-bit strobe/8-bit any edge';
end
if strcmpi(binaryValue(4), '1')
    h.PortReadType = '8-bit any bit falling edge';
end
if strcmpi(binaryValue(5), '1')
    h.PortReadType = '8-bit any bit rising edge';
end
if strcmpi(binaryValue(6), '1')
    h.PortControlsEvent = 1;
end
if strcmpi(binaryValue(7), '1')
    h.CapturePacketMethod = 'Logic input';
end
if strcmpi(binaryValue(8), '1')
    h.CapturePacketMethod = 'Character';
end
if strcmpi(binaryValue(9), '1')
    h.CapturePacketMethod = 'Word strob';
end
if strcmpi(binaryValue(10), '1')
    h.CapturePacketMethod = 'Any bit change';
end
if strcmpi(binaryValue(11), '1')
    h.NumberOfBits = 32;
end
if strcmpi(binaryValue(12), '1')
    h.NumberOfBits = 16;
end
if strcmpi(binaryValue(13), '1')
    h.NumberOfBits = 8;
end
if strcmpi(binaryValue(14), '1')
    h.NumberOfBits = 1;
end
if strcmpi(binaryValue(17), '1')
    h.BaudRate = 115200;
end
if strcmpi(binaryValue(18), '1')
    h.BaudRate = 57600;
end
if strcmpi(binaryValue(19), '1')
    h.BaudRate = 38400;
end
if strcmpi(binaryValue(20), '1')
    h.BaudRate = 19200;
end
if strcmpi(binaryValue(21), '1')
    h.BaudRate = 9600;
end
if strcmpi(binaryValue(22), '1')
    h.BaudRate = 2400;
end

function h = setFlagsCBAOUT(flags)

binaryValue = dec2bin(flags,22);

h.Audio = 0;
h.AnalogScale = 0;
h.ChanFunction = '';

if strcmpi(binaryValue(1), '1')
    h.ChanFunction = 'Waveform generator';
end
if strcmpi(binaryValue(2), '1')
    h.ChanFunction = 'Monitors spikes';
end
if strcmpi(binaryValue(3), '1')
    h.ChanFunction = 'Monitors SMP';
end
if strcmpi(binaryValue(4), '1')
    h.ChanFunction = 'Monitors LNC';
end
if strcmpi(binaryValue(5), '1')
    h.ChanFunction = 'Monitors raw';
end
if strcmpi(binaryValue(6), '1')
    h.AnalogScale = 'Waveform generator';
end
if strcmpi(binaryValue(7), '1')
    h.AnalogScale = 'Waveform generator';
end
if strcmpi(binaryValue(8), '1')
    h.AnalogScale = 'Waveform generator';
end
if strcmpi(binaryValue(9), '1')
    h.Audio = 1;
end

function h = setFlagsCBAINP(flags)

binaryValue = dec2bin(flags,22);

h.Audio = 0;
h.AnalogScale = 0;
h.ChanFunction = '';

if strcmpi(binaryValue(1), '1')
    h.ChanFunction = 'Waveform generator';
end
if strcmpi(binaryValue(2), '1')
    h.ChanFunction = 'Monitors spikes';
end
if strcmpi(binaryValue(3), '1')
    h.ChanFunction = 'Monitors SMP';
end
if strcmpi(binaryValue(4), '1')
    h.ChanFunction = 'Monitors LNC';
end
if strcmpi(binaryValue(5), '1')
    h.ChanFunction = 'Monitors raw';
end
if strcmpi(binaryValue(6), '1')
    h.AnalogScale = 'Waveform generator';
end
if strcmpi(binaryValue(7), '1')
    h.AnalogScale = 'Waveform generator';
end
if strcmpi(binaryValue(8), '1')
    h.AnalogScale = 'Waveform generator';
end
if strcmpi(binaryValue(9), '1')
    h.Audio = 1;
end


