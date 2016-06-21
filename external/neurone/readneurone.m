% readneurone()  -  Read data from a Mega NeurOne device and arrange it
%                   into a struct.
%
% Usage:    >> NEURONE = readneurone(dataPath, sessionPhaseNumber, chans)
%
% =======================================================================
% Inputs:
%   dataPath              - Path to directory containing NeurOne data 
%                           from a single measurement.
%
% Optional input arguments:
%   sessionPhaseNumber    - Defines which session phase will be read.
%                           Default value for sessionPhaseNumber is 1.
%   chans                 - Defines which channels will be read. The
%                           channel numbers have to be strings containing
%                           numeric values between 1 and the total number
%                           of channels. In addition, they are arranged in 
%                           vector form without brackets. For example, 
%                           if the total number of channels is 64, the 
%                           'chans' input argument can be '1,6,3,7' or 
%                           '1:12'. The default value for chans is an empty
%                           string which indicates that all channels will 
%                           be read. The channel numbers can be arranged 
%                           in a completely random order.
%
% =======================================================================
% Outputs:                           
%       NEURONE is a struct with the same fields as a standard EEG
%       structure under EEGLAB (see eeg_checkset.m).
%
% Used fields are:
%   NEURONE.filename
%   NEURONE.setname
%   NEURONE.srate
%   NEURONE.pnts
%   NEURONE.xmin
%   NEURONE.xmax
%   NEURONE.filepath
%   NEURONE.trials
%   NEURONE.nbchan
%   NEURONE.data
%   NEURONE.ref
%   NEURONE.times
%   NEURONE.comments
%   NEURONE.etc
%   NEURONE.subject
%   NEURONE.event (a struct containing all event data)
%   NEURONE.eventdescription
%   NEURONE.chanlocs (will be defined by EEGLAB if channel locations exist
%                for used channel names)
%   NEURONE.chaninfo
%
% The following fields were added as empty to ensure proper working:
%   NEURONE.icawinv
%   NEURONE.icaweights
%   NEURONE.icasphere
%   NEURONE.icaact
%
%
% ========================================================================
% NOTE:
% This file is part of the NeurOne data import plugin for EEGLAB.
% ========================================================================
% 
% Current version: 1.0.3.4 (2016-06-17)
% Author: Mega Electronics


function NEURONE = readneurone(dataPath, varargin)
%% Parse input arguments

p = inputParser;
p.addRequired('dataPath', @ischar);
p.addOptional('sessionPhaseNumber', 1 ,@isnumeric);
p.addOptional('chans','', @isstr);

p.parse(dataPath, varargin{:});
arglist = p.Results;

%%  Empty recording structure for output

NEURONE = {};

%% Read Session and Protocol xml-files

session = module_read_neurone_xml([dataPath 'Session.xml']);
protocol = module_read_neurone_xml([dataPath 'Protocol.xml']);

%% Obtain information from session and protocol

% Obtain sampling frequency (it is same for all channels)
srate = str2num(protocol.TableProtocol.SamplingFrequency);

% Get the total number of channels
nChannels = numel(protocol.TableInput);

NEURONE.comments = ['Protocol: ' session.TableSession.ProtocolName];
NEURONE.etc = strvcat( ['NeurOne Version: ' protocol.TableInfo.NeurOneVersion], ...
    ['Protocol File Revision: ' protocol.TableInfo.Revision], ...
    ['Session File Revision: ' session.TableInfo.Revision]);

NEURONE.subject = session.TablePerson.PersonID;

NEURONE.setname = session.TableSession.ProtocolName;

%% Check the session phase number

% Total number of sessions in recording
nSessionPhases = numel(session.TableSessionPhase);
disp(['Number of session phases: ' num2str(nSessionPhases)])
disp(['Current session phase number: ' num2str(arglist.sessionPhaseNumber)]);

% Check the validity of the sessionPhaseNumber taken from GUI
if arglist.sessionPhaseNumber>nSessionPhases
    disp('Current session phase number exceeds the total number of session phases.')
    disp('Defaulting to 1.')
    sessionPhaseNumber = 1;
elseif arglist.sessionPhaseNumber<=0
    disp('Invalid session phase number.')
    disp('Defaulting to 1.')
    sessionPhaseNumber = 1;
else
    sessionPhaseNumber = arglist.sessionPhaseNumber;
end

%% Get correct channel names and locations

% Set the channels for data acquisition
if isempty(arglist.chans)
    chans = 1:nChannels;
else
    chans = [str2num(arglist.chans)];
    chans = sort(chans); % ensure that the channels are in ascending order
end

% Check that no channel number exceeds the total number of channels
if any(chans>nChannels)
    error('One or more of the given channel numbers exceeds the total number of channels. See pop_readneurone help.')
elseif any(chans<1)
    error('One or more of the given channel numbers are invalid (zero or negative). See pop_readneurone help.')
end

% Get all channel names and corresponding input numbers 
inputNumbersAll = zeros(1,nChannels);
inputNumbersAll(1) = str2num(protocol.TableInput(1).InputNumber);
channelnames = protocol.TableInput(1).Name;

for n = 2:nChannels
    inputNumbersAll(n) = str2num(protocol.TableInput(n).InputNumber);
    channelnames = char(channelnames,protocol.TableInput(n).Name);
end

channelnames = cellstr(channelnames);

% The value of the highest input number
maxInput = max(inputNumbersAll);

channelnames = channelnames';

% Sort channel names in ascending order based on their input number
for k = 1:nChannels
    ii = 1;
    while ~(any(inputNumbersAll(k)==ii)) && ii<=max(inputNumbersAll)
          ii = ii+1;
    end
    sortIndex(ii) = ii; % Store the indices for future data arrangements
    channelNames(ii) = channelnames(k);
end

inputNumbersTrue = 1:maxInput;
inputNumbersTrue = inputNumbersTrue(find(sortIndex));

% Remove empty channel names:
% Find empty cells...
emptyCells = cellfun(@isempty,channelNames);
% ...and remove them.
channelNames(emptyCells) = [];

% Take only the names of the selected channels
channelNames = channelNames(chans);

fprintf('Looking up channel locations...\n')
chanlocs=struct('labels', channelNames');
% NEURONE.chanlocs = pop_chanedit(chanlocs);
NEURONE.chanlocs = chanlocs;

%% Obtain additional information about the dataset

measurementMode = protocol.TableInput(1).AlternatingCurrent;
deviceFilter = protocol.TableInput(1).Filter;
for n = 2:nChannels
    measurementMode = char(measurementMode,protocol.TableInput(n).AlternatingCurrent);
    deviceFilter = char(deviceFilter,protocol.TableInput(n).Filter);
end

measurementMode = cellstr(measurementMode);
deviceFilter = cellstr(deviceFilter);

for k = 1:nChannels
    if strcmp(measurementMode(k),'true')
        measurementMode(k) = {'AC'};
    else 
        measurementMode(k) = {'DC'};
    end
end

% Arrange measurementModes according to the input number
ii = 1;
for k = inputNumbersAll
    measurementMode(k) = measurementMode(ii);
    deviceFilter(k) = deviceFilter(ii);
    ii = ii+1;
end

measurementMode = measurementMode(inputNumbersTrue);
measurementMode = measurementMode(chans);
deviceFilter = deviceFilter(inputNumbersTrue);
deviceFilter = deviceFilter(chans);

measurementModeField = [char(channelNames(1)) ': ' char(measurementMode(1))];
deviceFilterField = [char(channelNames(1)) ': ' char(deviceFilter(1))];

for k = 2:length(chans)
    measurementModeField = strvcat(measurementModeField,[char(channelNames(k)) ': ' char(measurementMode(k))]);
    deviceFilterField = strvcat(deviceFilterField,[char(channelNames(k)) ': ' char(deviceFilter(k))]);
end

NEURONE.chaninfo.measurementMode = measurementModeField;
NEURONE.chaninfo.deviceFilter = deviceFilterField;

%% Preparing to read binary data
% If the size of the measurement in each session phase has exceeded the 
% file size, the data may have been split into several files named as 
% 1.bin, 2.bin etc. Usually there exists only one file: 1.bin. However, all
% this data needs to be read.

dataFiles = {}; % empty structure for data files

% Get all .bin files related to the chosen sessionPhaseNumber
sessionData = dir([dataPath num2str(sessionPhaseNumber) filesep '*.bin']);

ii = 1;
for k = 1:numel(sessionData); % total number of .bin files in current session phase data folder
    [path,filename,ext] = fileparts(sessionData(k).name);
    if ~isempty(regexpi(filename,'[123456789]'))
        dataFiles{ii,1} = [dataPath num2str(sessionPhaseNumber) filesep sessionData(k).name];
        ii = ii+1;
    end
end

%% Read NeurOne binary data and calibrate it

% Read data. The data of a specific channel is arranged in a row based on its input
% number.

data = readneuronedata(dataFiles, nChannels, chans);

% Preallocate memory to read maximum and minimum values for calibration
rawMinimum = zeros(1,maxInput);
rawMaximum = zeros(1,maxInput);
calibratedMinimum = zeros(1,maxInput);
calibratedMaximum = zeros(1,maxInput);

% Get all minimum and maximum values
for n = 1:nChannels
    rawMinimum(n) = str2num(protocol.TableInput(1,n).RangeMinimum);
    rawMaximum(n) = str2num(protocol.TableInput(1,n).RangeMaximum);
    calibratedMinimum(n) = str2num(protocol.TableInput(1,n).RangeAsCalibratedMinimum);
    calibratedMaximum(n) = str2num(protocol.TableInput(1,n).RangeAsCalibratedMaximum);
end

% Arrange them in the correct according to the input number
 ii = 1;
for k = inputNumbersAll
    rawMinimum(inputNumbersAll) = rawMinimum(ii);
    rawMaximum(inputNumbersAll) = rawMaximum(ii);
    calibratedMinimum(inputNumbersAll) = calibratedMinimum(ii);
    calibratedMaximum(inputNumbersAll) = calibratedMaximum(ii);
    ii = ii+1;
end

% Remove unused input numbers
rawMinimum = rawMinimum(inputNumbersTrue);
rawMaximum = rawMaximum(inputNumbersTrue);
calibratedMinimum = calibratedMinimum(inputNumbersTrue);
calibratedMaximum = calibratedMaximum(inputNumbersTrue);

% Calibrate channels
for n = 1:numel(chans)
    data(n,:) = calibratedMinimum(n) + ...
        (data(n,:)-rawMinimum(n)) / (rawMaximum(n)-rawMinimum(n)) * (calibratedMaximum(n)-calibratedMinimum(n));
end


%% Read events

fprintf('Loading events...\n');
event = readneuroneevents([dataPath num2str(sessionPhaseNumber) filesep]);
NEURONE.event = event;

% Writing event description
latency = strvcat('Latency:','The index number indicating the point in the data when the particular event has occurred');
type = strvcat('Type:','A descriptive name for the event. Will be generated depending on the type of the event as follows:',' ', ...
    'Event type', ...
    '----------', ...
    '1) Unknown', ... 
    '2) Stimulation', ...  
    '3) Video trigger', ...  
    '4) Mute trigger', ...    
    '5) 8-bit trigger', ...     
    '6) Comment', ' ', ...
    'Corresponding value of the type field', ...
    '-----------------------', ...
    '1) (Source port name) - Unknown', ...
    '2) (Source port name) - Stimulation', ...
    '3) (Source port name) - Video', ...
    '4) (Source port name) - Mute', ...
    '5) (8-bit trigger code, a value between 1-255)',...
    '6) (any user-entered comment related to the event)');

NEURONE.eventdescription = {latency type};

%% Store rest of the data
fprintf('Preparing output...\n');

% Basic dataset information:
NEURONE.srate = srate;
NEURONE.pnts = numel(data)/numel(chans);
NEURONE.xmin = 0;
NEURONE.xmax =(NEURONE.pnts-1)/srate;
NEURONE.trials = 1;
NEURONE.nbchan = length(chans);
NEURONE.data = data;
NEURONE.ref = 'common';
NEURONE.filepath = ''; % Will be generated by EEGLAB when saved
NEURONE.filename = ''; % Will be generated by EEGLAB when saved

% Additional dataset information:
NEURONE.icawinv = [];
NEURONE.icaweights = [];
NEURONE.icasphere = [];
NEURONE.icaact = [];

end