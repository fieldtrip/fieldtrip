function [sFile, ChannelMat] = in_fopen_nk(DataFile)
% IN_FOPEN_EDF: Open a Nihon Kohden file (.EEG / .PNT / .LOG / .21E)
%
% USAGE:  [sFile, ChannelMat] = in_fopen_nk(DataFile)

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2018 University of Southern California & McGill University
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
% Authors: Francois Tadel, 2017-2018
%          Inspired from NK2EDF, Teunis van Beelen, 2007-2017
%          and from the BIOSIG-toolbox http://biosig.sf.net/


%% ===== GET FILES =====
% Get base filename for all the files
[fPath, fBase] = bst_fileparts(DataFile);
BaseFile = bst_fullfile(fPath, fBase);
% EEG file (mandatory)
EegFile = [BaseFile '.eeg'];
if ~file_exist(EegFile)
    EegFile = [BaseFile '.EEG'];
    if ~file_exist(EegFile)
        error('Could not find .EEG file.');
    end
end
% PNT file (optional)
PntFile = [BaseFile '.pnt'];
if ~file_exist(PntFile)
    PntFile = [BaseFile '.PNT'];
    if ~file_exist(PntFile)
        disp('NK> Warning: Could not find .PNT file.');
        PntFile = [];
    end
end
% LOG file (optional)
LogFile = [BaseFile '.log'];
if ~file_exist(LogFile)
    LogFile = [BaseFile '.LOG'];
    if ~file_exist(LogFile)
        disp('NK> Warning: Could not find .LOG file.');
        LogFile = [];
    end
end
% 21E file (optional)
ElecFile = [BaseFile '.21e'];
if ~file_exist(ElecFile)
    ElecFile = [BaseFile '.21E'];
    if ~file_exist(ElecFile)
        disp('NK> Warning: Could not find .21E electrodes file.');
        LogFile = [];
    end
end


%% ===== READ EEG FILE =====
% Open file
fid = fopen(DataFile, 'rb');
if (fid == -1)
    error('Could not open EEG file.');
end
% Get deviceblock signature
hdr.device = fread(fid, [1 16], '*char');
hdr.version = get_header_version(hdr.device);
if (hdr.version == 0)
    error(['EEG deviceblock has unknown signature: "' hdr.device '"']);
end
% Get controlblock signature
fseek(fid, 129, 'bof');
hdr.control = fread(fid, [1 16], '*char');
if (get_header_version(hdr.control) == 0)
    error(['EEG controlblock has unknown signature: "' hdr.control '"']);
end
% Get waveformdatablock signature
fseek(fid, 6142, 'bof');
signature = fread(fid, [1 1], '*char');
if (signature ~= 1)
    error('waveformdatablock has wrong signature.');
end

% Get number of control blocks
fseek(fid, 145, 'bof');                        % Position:         0x0091
hdr.ctl_cnt = fread(fid, [1 1], 'uint8');
% Get a pointer to the extra block (only valid for newer file formats >= 1200 NK systems)
fseek(fid, 1006, 'bof');                       % Position:         0x03EE
hdr.ext_address = fread(fid, 1, 'uint32');     % extblock_address: VARIABLE
% Consider blocks of 100ms everywhere
hdr.record_duration = 0.1;

% Get all the pointers to all the blocks
for i = 1:hdr.ctl_cnt                               % Ctl block at:     0x0080
    % Get pointer to this block
    fseek(fid, 146 + (i-1) * 20, 'bof');            % Position:         0x0092   (1st ctl block)
    hdr.ctl(i).address = fread(fid, 1, 'uint32');   % ctlblock_address: 0x0400
    % Get number of data blocks   
    fseek(fid, hdr.ctl(i).address + 23, 'bof');     % Position:         0x0417
    hdr.ctl(i).data_cnt = fread(fid, [1 1], 'uint8');
    hdr.ctl(i).data = repmat(struct(), 0);
    % Loop on data blocks
    for j = 1:hdr.ctl(i).data_cnt
        % Read data address
        fseek(fid, hdr.ctl(i).address + ((j-1) * 20) + 18, 'bof');            % Position:       0x0412   (1st ctl block)
        dataAddr = fread(fid, 1, 'uint32');                                   % Data block at:  0x17FE
        % Add data block only if it points to a valid address in the file
        if (dataAddr > 0)            
            % Save address in structure
            id = length(hdr.ctl(i).data) + 1;
            hdr.ctl(i).data(id).address = dataAddr;
            % Read block start timestamp
            fseek(fid, dataAddr + 5, 'bof');
            timeH = str2double(fread(fid, [1 2], '*char'));
            timeM = str2double(fread(fid, [1 2], '*char'));
            timeS = str2double(fread(fid, [1 2], '*char'));
            hdr.ctl(i).data(id).timestamp = 60*60*timeH + 60*timeM + timeS;
            % Read the sampling rate
            fseek(fid, dataAddr + 26, 'bof');              % Position:   0x1818
            hdr.ctl(i).data(id).sample_rate = bitand(fread(fid, 1, 'uint16'), hex2dec('3fff'));    % Nihon-Kohden int16 format

            % Get the channel order from here (older system)
            switch (hdr.version)
                case 1    % Older NK systems: 1100, 2100
                    % Read the block information: number of records
                    fseek(fid, dataAddr + 28, 'bof');      % Position:   0x181A
                    hdr.ctl(i).data(id).num_records = fread(fid, 1, 'uint32');
                    % Read number of channels
                    fseek(fid, dataAddr + 38, 'bof');      % Position:   0x1824
                    hdr.ctl(i).data(id).num_channels = fread(fid, 1, 'uint8') + 1;  % +1 for the STIM channel
                    % Read channel order
                    for iChan = 1:(hdr.ctl(i).data(id).num_channels - 1)    % -1 because the STIM channel is not listed here
                        fseek(fid, dataAddr + 39 + (iChan - 1) * 10, 'bof');  % Position:   0x1825, 0x182F, 0x1839, ...
                        hdr.ctl(i).data(id).channel_list(iChan) = fread(fid, 1, 'uint8') + 1;
                    end
                    % Define pointer to the beginning of the recordings
                    hdr.ctl(i).data(id).rec_address = hdr.ctl(i).data(id).address + 39 + (hdr.ctl(i).data(id).num_channels - 1) * 10;   % -1 because the STIM channel is not listed here
                    % Compute number of samples
                    hdr.ctl(i).data(id).num_samples = hdr.ctl(i).data(id).num_records * hdr.ctl(i).data(id).sample_rate * hdr.record_duration;
                case 2    % Newers NK systems: 1200
                    % Channel order read in the extended blocks
                    hdr.ctl(i).data(id).num_records = [];    % TODO: DON'T KNOW WHERE TO GET THIS FROM IN THE FILE
                    hdr.ctl(i).data(id).num_samples = [];    % TODO: DON'T KNOW HOW TO COMPUTE THIS (NOW GUESSING IT FROM FILE SIZE)
            end
        end
    end
end

% Read information from additional blocks (newer systems)
switch (hdr.version)
    case 1    % Older NK systems: 1100, 2100
        % Not needed here, corresponding pointer is 0x0000 in the file
    case 2    % Newers NK systems: 1200
        % Only supported for one control block + one data block
        if (length(hdr.ctl) > 1) || (length(hdr.ctl(1).data) > 1)
            error(['This reader supports only recordings file with one data segment.' 10 ...
                   'If you are interested in reading files with multiple data segments,' 10 ...
                   'please contact us through the Brainstorm user forum.']);
        end
        
        % TODO: Probably needs a loop on control and data blocks here => Need example files
        i = 1;
        id = 1;
        
        % Reading the extended block address (2nd pointer)
        % (hdr.ext_address + 17) = UINT8 Number of blocks ?
        fseek(fid, hdr.ext_address + 18, 'bof');
        hdr.ctl(i).extblock2_address = fread(fid, 1, 'uint32');

        % Reading the extended block address (3rd pointer)
        % (hdr.ext_address2 + 17) = UINT8 Number of blocks ?
        fseek(fid, hdr.ctl(i).extblock2_address + 20, 'bof');
        hdr.ctl(i).data(id).extblock3_address = fread(fid, 1, 'uint32');
        
% Suggestion from V. Gnatkovsky: not working with some of the files...
%         % Read the block information: number of records 
%         fseek(fid, hdr.ctl(i).data(id).extblock3_address + 44, 'bof');
%         hdr.ctl(i).data(id).num_records = fread(fid, 1, 'uint32');
%         % Compute number of samples
%         hdr.ctl(i).data(id).num_samples = hdr.ctl(i).data(id).num_records * hdr.ctl(i).data(id).sample_rate * hdr.record_duration;

        % Reading number of channels
        fseek(fid, hdr.ctl(i).data(id).extblock3_address + 68, 'bof');
        hdr.ctl(i).data(id).num_channels = fread(fid, 1, 'uint16') + 1;   % +1 for the STIM channel
        % Read channel order
        for iChan = 1:(hdr.ctl(i).data(id).num_channels - 1)  % -1 because the STIM channel is not listed here
            fseek(fid, hdr.ctl(i).data(id).extblock3_address + 72 + (iChan-1) * 10, 'bof');
            hdr.ctl(i).data(id).channel_list(iChan) = fread(fid, 1, 'uint16') + 1;
        end
        % Define pointer to the beginning of the recordings
        hdr.ctl(i).data(id).rec_address = hdr.ctl(i).data(id).extblock3_address + 72 + (hdr.ctl(i).data(id).num_channels - 1) * 10;   % -1 because the STIM channel is not in this list
end

% Get last position in the file
fseek(fid, 0, 'eof');
lastpos = ftell(fid);
% Close file
fclose(fid);


%% ===== GUESS FILE PROPERTIES =====
% TODO: Current limitation: allow only files with single control blocks
if (length(hdr.ctl) ~= 1)
    error(['Files with more than one control block are currently not supported.' 10 ...
           'Please post a message on the Brainstorm forum if you need this feature to be enabled.']);
end
% TODO: Current limitation: Multiple data blocks are not supported yet with the new file format (just need an example dataset)
if (length(hdr.ctl(1).data) > 1) && (hdr.version == 2)
    error(['Newer files with more than one data block are currently not supported yet (systems NK EEG-1200A V01.00).' 10 ...
           'Please post a message on the Brainstorm forum if you need this feature to be enabled.']);
end
% TODO: Current limitation: Multiple data blocks must have the same properties
if (length(hdr.ctl(1).data) > 1) && (...
        any([hdr.ctl(1).data.num_channels] ~= hdr.ctl(1).data(1).num_channels) || ...
        any([hdr.ctl(1).data.sample_rate] ~= hdr.ctl(1).data(1).sample_rate))
    error('Files with more than one data block must have the same number of channels and the same sampling rate.');
end
% Copy shared fields to the central header
hdr.sample_rate  = hdr.ctl(1).data(1).sample_rate;
hdr.num_channels = hdr.ctl(1).data(1).num_channels;


%% ===== READ LOG FILE =====
if ~isempty(LogFile)
    % Open file
    fid = fopen(LogFile, 'rb');
    if (fid == -1)
        error('Could not open LOG file');
    end
    % Get file signature
    device = fread(fid, [1 16], '*char');
    if (get_header_version(device) == 0)
        error(['LOG file has unknown signature: "' device '"']);
    end
    % Get log blocks 
    fseek(fid, 145, 'bof');
    n_logblocks = fread(fid, 1, 'uint8');
    % Initializations
    total_logs = 0;
    
    % Loop on log blocks
    for i = 1:n_logblocks
        % Read number of logs in this block
        fseek(fid, 146 + ((i-1) * 20) , 'bof');
        logblock_address = fread(fid, 1, 'uint32');
        fseek(fid, logblock_address + 18, 'bof');
        n_logs = fread(fid, 1, 'uint8');
        % Initialization
        fseek(fid, logblock_address + 20, 'bof');
        hdr.logs(i).label = cell(1, n_logs);
        hdr.logs(i).time  = zeros(1, n_logs);
        % Read all the events
        for j = 1:n_logs
            hdr.logs(i).label{j} = str_clean(fread(fid, [1 20], '*char'));
            timeH = str2double(fread(fid, [1 2], '*char'));
            timeM = str2double(fread(fid, [1 2], '*char'));
            timeS = str2double(fread(fid, [1 2], '*char'));
            hdr.logs(i).time(j) = 60*60*timeH + 60*timeM + timeS;
            hdr.logs(i).label2{j} = fread(fid, [1 19], '*char');
            % Compute time stamp
            timeH = str2double(hdr.logs(i).label2{j}(8:9));
            timeM = str2double(hdr.logs(i).label2{j}(10:11));
            timeS = str2double(hdr.logs(i).label2{j}(12:13));
            hdr.logs(i).timestamp(j) = 60*60*timeH + 60*timeM + timeS;
        end
            
        % Read sub-events
        try
            % Read number of sub-logs
            fseek(fid, 146 + (((i-1) + 22) * 20) , 'bof');
            sublogblock_address = fread(fid, 1, 'uint32');
            fseek(fid, sublogblock_address + 18, 'bof');
            n_sublogs = fread(fid, 1, 'uint8');
            % Read sub-logs
            if (n_sublogs == n_logs)
                fseek(fid, sublogblock_address + 20, 'bof');
                for j = 1:n_logs
                    hdr.logs(i).sublog{j} = fread(fid, [1 45], '*char');
                    hdr.logs(i).time(j) = hdr.logs(i).time(j) + str2double(['0.' hdr.logs(i).sublog{j}(25:30)]);
                end
            end
        catch
            disp('NK> Could not read sub-events.');
        end
        total_logs = total_logs + n_logs;
    end
    % Close file
    fclose(fid);
else
    hdr.logs = [];
end


%% ===== CHANNELS: READ 21E FILE =====
% Read the channel names 
ChannelMat = in_channel_nk(ElecFile, hdr.version);

% Gains are fixed for the list of channels: 
%   - uV for channels: 1-42, 75, 76, 79-256, 257:1096 (new systems)
%   - mV for all the others
%   - Calibration = (Physical_max - Physical_min) ./ (Digital_max - Digital_min)
iChanMicro = [1:42, 75, 76, 79:256, 257:1096];
chanGains = 1e-3 * ones(1,1096) * ((12002.56+12002.9) / (32767 + 32768));
chanGains(iChanMicro) = 1e-6 * ((3199.902+3200) / (32767 + 32768));

% Keep only the channels saved in the file
iSelChannels = hdr.ctl(1).data(id).channel_list;
ChannelMat.Channel = ChannelMat.Channel(iSelChannels);
hdr.channel_gains  = chanGains(iSelChannels);

% Add last channel: markers/events
iChanStim = length(ChannelMat.Channel) + 1;
ChannelMat.Channel(iChanStim).Name    = 'Events';
ChannelMat.Channel(iChanStim).Type    = 'STIM';
hdr.channel_gains(iChanStim) = 1;

% Digital offset for the channel calibration (only the last one is a digital channel => no offset)
hdr.channel_digoffset = [32768 .* ones(1, length(ChannelMat.Channel)), 0];

% Make sure the names of the channels are unique
for i = length(ChannelMat.Channel):-1:2
    if ismember(ChannelMat.Channel(i).Name, {ChannelMat.Channel(1:i-1).Name})
        ChannelMat.Channel(i).Name = [ChannelMat.Channel(i).Name, '_', num2str(i)];
        ChannelMat.Channel(i).Type = 'MISC';
    end
end
% % Detect what is EEG or SEEG
% iTypeEeg = find(strcmpi({ChannelMat.Channel.Type}, 'EEG'));
% iNameEeg = find(ismember({ChannelMat.Channel.Name}, {'FP1', 'FP2', 'CZ', 'FCZ', 'PZ', 'OZ'}));
% if ~ismember({ChannelMat.Channel.Name}, 'F1')
%     iNameEeg = [iNameEeg, find(ismember({ChannelMat.Channel.Name}, {'F3', 'F4', 'F7', 'F8'}))];
% end
% if ~ismember({ChannelMat.Channel.Name}, 'C1')
%     iNameEeg = [iNameEeg, find(ismember({ChannelMat.Channel.Name}, {'C3', 'C4'}))];
% end
% if ~ismember({ChannelMat.Channel.Name}, 'P1')
%     iNameEeg = [iNameEeg, find(ismember({ChannelMat.Channel.Name}, {'P3', 'P4'}))];
% end
% if ~ismember({ChannelMat.Channel.Name}, 'O5')
%     iNameEeg = [iNameEeg, find(ismember({ChannelMat.Channel.Name}, {'O1', 'O2', 'O3', 'O4'}))];
% end
% if ~ismember({ChannelMat.Channel.Name}, 'T1')
%     iNameEeg = [iNameEeg, find(ismember({ChannelMat.Channel.Name}, {'T3', 'T4', 'T5', 'T6'}))];
% end
% % Mark all the other channels as SEEG
% iSeeg = setdiff(iTypeEeg, iNameEeg);
% [ChannelMat.Channel(iSeeg).Type] = deal('SEEG');


%% ===== READ PNT FILE =====
if ~isempty(PntFile)
    % Open file
    fid = fopen(PntFile, 'rb');
    if (fid == -1)
        error('Could not open LOG file');
    end
    % Get file signature
    device = fread(fid, [1 16], '*char');
    if (get_header_version(device) == 0)
        error(['PNT file has unknown signature: "' device '"']);
    end
    % Read patient info: Id
    fseek(fid, 1540, 'bof');
    hdr.patient.Id = str_clean(fread(fid, [1 10], '*char'));
    % Read patient info: Name
    fseek(fid, 1582, 'bof');
    hdr.patient.Name = str_clean(fread(fid, [1 20], '*char'));
    % Read patient info: Sex
    fseek(fid, 1610, 'bof');
    hdr.patient.Sex = str_clean(fread(fid, [1 6], '*char'));
    % Read patient info: Birthday
    fseek(fid, 1632, 'bof');
    hdr.patient.Birthday = fread(fid, [1 10], '*char');
    % Read recordings date
    fseek(fid, 64, 'bof');
    numDate = sscanf(fread(fid, [1 14], '*char'), '%04u%02u%02u%02u%02u%02u');
    hdr.startdate = sprintf('%02d/%02d/%04d', numDate(3), numDate(2), numDate(1));
    hdr.starttime = sprintf('%02d:%02d:%02d', numDate(4), numDate(5), numDate(6));
    % Close file
    fclose(fid);
end


%% ===== CREATE BRAINSTORM SFILE STRUCTURE =====
% Initialize returned file structure
sFile = db_template('sfile');
% Add information read from header
sFile.byteorder = 'n';
sFile.filename  = DataFile;
sFile.format    = 'EEG-NK';
sFile.device    = ['Nihon Kohden ' hdr.device];
sFile.comment   = fBase;
sFile.condition = [];
% Epochs
nEpochs = length(hdr.ctl(1).data);
sFile.epochs = repmat(db_template('epoch'), 1, nEpochs);
for i = 1:nEpochs
    % If the number of samples is not known (new format): guess it from the file or block size
    if isempty(hdr.ctl(1).data(i).num_samples)
        % There is a block after: use the address of the next block
        if (i < nEpochs)
            end_address = hdr.ctl(1).data(i+1).address - 1;
        % Last block: use the file size
        else
            end_address = lastpos;
        end
        % Compute from the 
        hdr.ctl(1).data(i).num_samples = floor((end_address - hdr.ctl(1).data(i).rec_address) / hdr.num_channels / 2);  % /2 because we are counting int16 values
    end
    sFile.epochs(i).samples = [0, hdr.ctl(1).data(i).num_samples - 1];
    sFile.epochs(i).times   = sFile.epochs(i).samples ./ hdr.sample_rate;
    sFile.epochs(i).label   = sprintf('Block #%d', i);
    sFile.epochs(i).nAvg    = 1;
    sFile.epochs(i).select  = 1;
    sFile.epochs(i).bad     = 0;
    % Cumulated time of the begginning of the epoch
    epochLength(i) = hdr.ctl(1).data(i).num_samples ./ hdr.sample_rate;
end
cumTime = [0, cumsum(epochLength(1:end-1))];
% Consider that the sampling rate of the file is the sampling rate of the first signal
sFile.prop.sfreq   = hdr.sample_rate;
sFile.prop.samples = [min([sFile.epochs.samples]), max([sFile.epochs.samples])];
sFile.prop.times   = [min([sFile.epochs.times]),   max([sFile.epochs.times])];
sFile.prop.nAvg    = 1;
% No info on bad channels
sFile.channelflag = ones(hdr.num_channels,1);
% Save full header in the file link
sFile.header = hdr;
% Acquisition date
sFile.acq_date = str_date(hdr.startdate);


%% ===== EVENTS =====
if ~isempty(hdr.logs)
    % Get all the event types
    evtList = hdr.logs(1).label;
    % Events list
    [uniqueEvt, iUnique] = unique(evtList);
    uniqueEvt = evtList(sort(iUnique));
    % Initialize events list
    sFile.events = repmat(db_template('event'), 1, length(uniqueEvt));
    % Build events list
    for iEvt = 1:length(uniqueEvt)
        % Find all the occurrences of this event
        iOcc = find(strcmpi(uniqueEvt{iEvt}, evtList));
        % Concatenate all times
        t = hdr.logs(1).time(iOcc);
        % Detect epoch
        sFile.events(iEvt).epochs = 1 + 0*t(1,:);
        for k = 1:length(iOcc)
            % iEpoch = find(hdr.logs(1).timestamp(iOcc(k)) <= [hdr.ctl(1).data.timestamp, Inf], 1) - 1;
            iEpoch = find(t(k) >= cumTime, 1, 'last');
            if (iEpoch > 1)
                sFile.events(iEvt).epochs(k) = iEpoch;
                t(k) = t(k) - cumTime(iEpoch);
            end
        end
        % Set event
        sFile.events(iEvt).label   = str_clean(uniqueEvt{iEvt});
        sFile.events(iEvt).select  = 1;
        sFile.events(iEvt).times   = t;
        sFile.events(iEvt).samples = round(t .* sFile.prop.sfreq);
    end
end


end



%% ===== CHECK DEVICE =====
function ver = get_header_version(str)
    % Older NK systems
    if ismember(str, {...
        'EEG-1100A V01.00', ...
        'EEG-1100B V01.00', ...
        'EEG-1100C V01.00', ...
        'QI-403A   V01.00', ...
        'QI-403A   V02.00', ...
        'EEG-2100  V01.00', ...
        'EEG-2100  V02.00', ...
        'DAE-2100D V01.30', ...
        'DAE-2100D V02.00', ...
        'EEG-1100A V02.00', ...
        'EEG-1100B V02.00', ...
        'EEG-1100C V02.00'})
        ver = 1;
        
    % Newer NK systems (>= 2015)
    elseif ismember(str, {...
        'EEG-1200A V01.00'})
        ver = 2;
        
    else
        ver = 0;
    end
end

%% ===== CLEAN STRINGS =====
function s = str_clean(s)
    % Stop string at first termination
    iNull = find(s == 0, 1);
    if ~isempty(iNull)
        s(iNull:end) = [];
    end
    % Remove weird characters
    s(~ismember(s, '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-,:;.*+=?!<>''"`&%$()[]{}/\_@ ·¡‡¿‚¬‰ƒ„√Â≈Ê∆Á«È…Ë»Í ÎÀÌÕÏÃÓŒÔœÒ—Û”Ú“Ù‘ˆ÷ı’¯ÿúåﬂ˙⁄˘Ÿ˚€¸‹')) = [];
    % Remove useless spaces
    s = strtrim(s);
end
    

