% readneuroneevents()   -   Read events from a Mega NeurOne device.
%
% Usage:    >> event = readneuroneevents(dataPath)
%
% =======================================================================
% Input:
%       dataPath        - Direct path for the folder containing the event
%                         data (file 'events.bin') related to the current 
%                         session phase number.
%
% Output:
%       event           - A struct created according to the event structure
%                         under EEGLAB. The used fields are type and
%                         latency. The field 'type' will be generated
%                         depending on the type of the event as follows:
%
%          Event type               Value of the type field
%          ----------               -----------------------
%          Unknown                 (Source port name) - Unknown
%          Stimulation             (Source port name) - Stimulation
%          Video trigger           (Source port name) - Video
%          Mute trigger            (Source port name) - Mute
%          8-bit trigger           (8-bit trigger code, a value between
%                                   1-255)
%          Comment                 (any user-entered comment related to
%                                   the event)
%          ClockSourceChange       ClockSourceChange (when the source of 
%                                  SyncBox clock is changed)
%
% =======================================================================
% Additional information:
%   The size of one event structure 'events.bin' file is 88 bytes. If user-
%   entered comments are used, they will be read from file 'eventData.bin' 
%   which holds comment for each event in Unicode (2 bytes/character). See
%   'NeurOne Data Format' documentation for more info about the event
%   structure.
%
% ========================================================================
% NOTE:
% This file is part of the NeurOne data import plugin for EEGLAB.
% ========================================================================
% 
% Current version: 1.0.3.4 (2016-06-17)
% Author: Mega Electronics


function event = readneuroneevents(dataPath)

% Get the total number of events
eventsTmp = dir([dataPath filesep 'events.bin']);
nEvents = eventsTmp.bytes/88;

event = {}; % empty structure for event data

% Read events.bin (see NeurOne Data Format doc)
events = fopen([dataPath filesep 'events.bin'],'rb');
for k = 1:nEvents
    % Read the whole event structure
    Revision = fread(events,1,'int32');          
    RFU = fread(events,1,'int32');
    Type = fread(events,1,'int32');  
    SourcePort = fread(events,1,'int32'); 
    ChannelNumber = fread(events,1,'int32');
    Code = fread(events,1,'int32');
    StartSampleIndex = fread(events,1,'uint64');
    StopSampleIndex = fread(events,1,'uint64');
    DescriptionLength = fread(events,1,'uint64');
    DescriptionOffset = fread(events,1,'uint64');
    DataLength = fread(events,1,'uint64');
    DataOffset = fread(events,1,'uint64');
    TimeStamp = fread(events,1,'double');
    MainUnitIndex = fread(events,1,'int32');
    RFU = fread(events,1,'int32');
    
    % Check if the data format has changed
    if Revision > 6
        warning(strcat('This reader does not support the revision of events.bin (', ...
            num2str(Revision), '). Please contact mega@megaemg.com for an update.'))
    end

    % Determine the source port
    switch SourcePort
        case 0
            SourcePort = 'N/A';
        case 1
            SourcePort = 'A';
        case 2
            SourcePort = 'B';
        case 3
            SourcePort = 'EightBit';
		case 4
            SourcePort = 'Syncbox Button';
        case 5
            SourcePort = 'SyncBox EXT';
        case 6
            SourcePort = 'Software';
        otherwise
            SourcePort = 'Unknown';
    end
    
    % Determine the type of the event
    switch Type
        case 0
            Type = [SourcePort ' - N/A'];
        case 1
            Type = [SourcePort ' - Stimulation'];
        case 2
            Type = [SourcePort ' - Video'];
        case 3
            Type = [SourcePort ' - Mute'];
        case 4
            Type = num2str(Code);
        case 5
            Type = [SourcePort ' - Out' ];
        case 6
            % User-entered comments will be read from file eventData.bin
            fid = fopen([dataPath filesep 'eventData.bin'],'rb');
            offset = DataOffset/2;
            length = DataLength/2;
            comments = fread(fid,[1 offset],'int16');
            comments = fread(fid,[1 length],'int16');
            Type = char(comments);
            fclose(fid);
        case 1001
            Type = [SourcePort ' - ClockSourceChange'];
		otherwise
            Type = 'Unknown';
    end
     
    % Store the obtained data
    event(1,k).latency = StartSampleIndex;
    event(1,k).type = Type;
   
end
fclose(events);

end
