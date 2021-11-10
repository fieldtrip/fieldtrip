function [sFile, ChannelMat] = in_fopen_manscan(DataFile)
% IN_FOPEN_MANSCAN: Open a MANSCAN file (continuous recordings)
%
% USAGE:  [sFile, ChannelMat] = in_fopen_manscan(DataFile)

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2012-2018
        

%% ===== READ HEADER =====
% Get text file (.mbi)
MbiFile = strrep(DataFile, '.mb2', '.mbi');
% If doesn't exist: error
if ~file_exist(MbiFile)
    error('Cannot open file: missing text file .mbi');
end
% Initialize header
iEpoch = 1;
hdr.epoch(iEpoch).Comment = {};
hdr.epoch(iEpoch).Channel = [];
hdr.Events = [];
curBlock = '';
% Read file line by line
fid = fopen(MbiFile,'r');
while(1)
    % Reached the end of the file: exit the loop
    if feof(fid)
        break; 
    end
    % Read one line
    buf = fgetl(fid);
    if isempty(buf)
        curBlock = '';
        continue;
    end
    % Split line based on space characters
    splitBuf = str_split(buf, [0 9 32]);
    if isempty(splitBuf)
        continue;
    end
    % Current block
    if ~isempty(curBlock)
        switch curBlock
            case 'event'
                % Ignore MANSCAN auto events
                if any(strcmpi(evtProp, 'Channel'))
                    continue;
                end
                % New event
                iEvt = length(hdr.Events) + 1;
                % Read event name and epoch
                hdr.Events(iEvt).Name   = splitBuf{1};
                hdr.Events(iEvt).iEpoch = iEpoch;
                % Read each property
                for iProp = 2:length(evtProp)
                    if strcmpi(evtProp{iProp}, 'Channel')
                        hdr.Events(iEvt).Name = [hdr.Events(iEvt).Name '_' splitBuf{iProp}];
                    elseif strcmpi(evtType{iProp}, 'String')
                        hdr.Events(iEvt).(evtProp{iProp}) = splitBuf{iProp};
                    else
                        hdr.Events(iEvt).(evtProp{iProp}) = str2num(splitBuf{iProp});
                    end
                end
                
            case 'channel'
                % New channel
                iChan = length(hdr.epoch(iEpoch).Channel) + 1;
                % Read channel name
                hdr.epoch(iEpoch).Channel(iChan).Name = splitBuf{1};
                % Read each property
                for iProp = 2:length(chanProp)
                    if strcmpi(chanProp{iProp}, 'UnitDisplayToFile')
                        propValue = str2num(splitBuf{iProp});
                    else
                        propValue = splitBuf{iProp};
                    end
                    hdr.epoch(iEpoch).Channel(iChan).(chanProp{iProp}) = propValue;
                end
        end
        % Next line
        continue;
    end
    % Detect keyword
    switch lower(splitBuf{1})
        case 'datatypeid'
            iEpoch = iEpoch + 1;
            hdr.epoch(iEpoch).Comment = {};
            hdr.epoch(iEpoch).Channel = [];
        case 'comment'
            hdr.epoch(iEpoch).Comment{end+1} = strtrim(strrep(buf, splitBuf{1}, ''));
        case 'event'
            curBlock = 'event';
            evtProp = splitBuf;
            evtType = str_split(fgetl(fid), [0 9 32]);
        case 'channel'
            curBlock = 'channel';
            chanProp = splitBuf;
        case 'offsetdisplayexperiment'
            hdr.epoch(iEpoch).OffsetDisplay = str2num(splitBuf{2});
        case 'experimenttime'
            hdr.epoch(iEpoch).ExperimentTime = strtrim(strrep(buf, splitBuf{1}, ''));
        case 'worddatafile'
            % Read channel order
            hdr.epoch(iEpoch).ChannelOrder = str_split(fgetl(fid), [0 9 32]);
            % Read record info
            datainfo = str_split(fgetl(fid), [0 9 32]);
            hdr.epoch(iEpoch).StartData1 = str2num(datainfo{1});
            hdr.epoch(iEpoch).StartData2 = str2num(datainfo{2});
            hdr.epoch(iEpoch).Sweeps     = str2num(datainfo{3});
            hdr.epoch(iEpoch).BinFile    = datainfo{4};
    end
end


%% ===== CREATE BRAINSTORM SFILE STRUCTURE =====
% Initialize returned file structure
sFile = db_template('sfile');
% Add information read from header
sFile.byteorder  = 'l';
sFile.filename   = DataFile;
sFile.format = 'EEG-MANSCAN';
sFile.device = 'MANSCAN';
% Comment: short filename
[tmp__, sFile.comment, tmp__] = bst_fileparts(DataFile);
% Consider that the sampling rate of the file is the sampling rate of the first signal
sFile.prop.sfreq = 256;
sFile.prop.nAvg  = 1;
% No info on bad channels
sFile.channelflag = ones(length(hdr.epoch(iEpoch).ChannelOrder), 1);
% Acquisition date
try
    sFile.acq_date = str_date(hdr.epoch.ExperimentTime, 'mm/dd/yy');
catch
end

%% ===== EPOCHS =====
if (length(hdr.epoch) <= 1)
    sFile.prop.times   = [0, hdr.epoch(1).Sweeps-1] ./ sFile.prop.sfreq;
else
    % Build epochs structure
    for iEpoch = 1:length(hdr.epoch)
        sFile.epochs(iEpoch).label   = sprintf('Epoch #%02d', iEpoch);
        sFile.epochs(iEpoch).times   = [0, hdr.epoch(iEpoch).Sweeps-1] ./ sFile.prop.sfreq;
        sFile.epochs(iEpoch).nAvg    = 1;
        sFile.epochs(iEpoch).select  = 1;
        sFile.epochs(iEpoch).bad     = 0;
        sFile.epochs(iEpoch).channelflag = [];
        
        % Check if all the epochs have the same channel list
        if (iEpoch > 1) && ~isequal(hdr.epoch(iEpoch).ChannelOrder, hdr.epoch(1).ChannelOrder)
            error('Channel list must remain constant across epochs.');
        end
    end
    % Extract global min/max for time and samples indices
    sFile.prop.times   = [min([sFile.epochs.times]),   max([sFile.epochs.times])];
end


%% ===== CREATE EMPTY CHANNEL FILE =====
nChannels = length(hdr.epoch(1).ChannelOrder);
ChannelMat = db_template('channelmat');
ChannelMat.Comment = [sFile.device ' channels'];
ChannelMat.Channel = repmat(db_template('channeldesc'), [1, nChannels]);
hdr.Gains = [];
% For each channel
for iChan = 1:nChannels
    % Find channel in description list
    chName = hdr.epoch(1).ChannelOrder{iChan};
    iDesc = find(strcmpi({hdr.epoch(1).Channel.Name}, chName));
    sDesc = hdr.epoch(1).Channel(iDesc);
    % Type
    if isfield(sDesc, 'IsEEG') && ~isempty(sDesc.IsEEG)
        if strcmpi(sDesc.IsEEG, 'Yes')
            ChannelMat.Channel(iChan).Type = 'EEG';
        else
            ChannelMat.Channel(iChan).Type = 'Misc';
        end
    elseif isfield(sDesc, 'IsEEG')
        ChannelMat.Channel(iChan).Type = 'EEG REF';
    else
        ChannelMat.Channel(iChan).Type = 'EEG';
    end
    % Create structure
    ChannelMat.Channel(iChan).Name    = chName;
    ChannelMat.Channel(iChan).Loc     = [0; 0; 0];
    ChannelMat.Channel(iChan).Orient  = [];
    ChannelMat.Channel(iChan).Weight  = 1;
    ChannelMat.Channel(iChan).Comment = '';
    % Check that this channel has a gain
    if isempty(hdr.epoch(1).Channel(iDesc).UnitDisplayToFile)
        hdr.Gains(iChan,1) = 1;
    else
        hdr.Gains(iChan,1) = hdr.epoch(1).Channel(iDesc).UnitDisplayToFile;
    end
end


%% ===== FORMAT EVENTS =====
if ~isempty(hdr.Events)
    % Get all the epochs names
    uniqueEvt = unique({hdr.Events.Name});
    % Create one category for each event
    for iEvt = 1:length(uniqueEvt)
        % Get all the occurrences
        iOcc = find(strcmpi({hdr.Events.Name}, uniqueEvt{iEvt}));
        % Get the samples for all the occurrences
        sample   = [hdr.Events(iOcc).BeginSample];
        duration = [hdr.Events(iOcc).DurationSamples];
        % Extended event = duration is not 1 for all the markers
        if ~all(duration == 1)
            sample = [sample; sample + duration];
        end
        % Create event structure
        sFile.events(iEvt).label    = uniqueEvt{iEvt};
        sFile.events(iEvt).times    = sample ./ sFile.prop.sfreq;
        sFile.events(iEvt).epochs   = [hdr.Events(iOcc).iEpoch];
        sFile.events(iEvt).select   = 1;
        sFile.events(iEvt).channels = cell(1, size(sFile.events(iEvt).times, 2));
        sFile.events(iEvt).notes    = cell(1, size(sFile.events(iEvt).times, 2));
    end
end
% Save file header
sFile.header = hdr;
    
    

