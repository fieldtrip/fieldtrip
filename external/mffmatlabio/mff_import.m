% mff_import - import MFF file to EEGLAB structure. This function calls
%                 all other function to import MFF xml files, including
%                 events, channels and channel coordinates.
%
% Usage:
%   EEG = mff_import(mffFile);
%
% Input:
%  mffFile - filename/foldername for the MFF file (MFF file/folder must
%            already exist)
% Output:
%  EEG      - EEGLAB structure
%  data     - channels x sample data
%  events   - events
%  chanlocs - channel locations
%  mff      - misceleneous event information
%
% Note: This function imports the "code" MFF event information into the
%       EEGLAB type field. If you want to use other MFF event information
%       extract epoch, use the pop_importmff function.

% This is the list of known types
% kMFF_RT_Unknown
% kMFF_RT_Any
% kMFF_RT_MFFFile
% kMFF_RT_Signal
% kMFF_RT_EventTrack
% kMFF_RT_Epochs
% kMFF_RT_Subject
% kMFF_RT_History
% kMFF_RT_Info
% kMFF_RT_InfoN
% kMFF_RT_Categories
% kMFF_RT_JTFCategories
% kMFF_RT_SensorLayout
% kMFF_RT_Coordinates
% kMFF_RT_Photogrammetry^
% kMFF_RT_PNSSet
% kMFF_RT_MovieSyncs
% kMFF_RT_Fields
% kMFF_RT_Notes
% kMFF_RT_Montage
% kMFF_RT_DipoleSet
% kMFF_RT_PhoticStim
% kMFF_RT_GTENModulationConfiguratonFile
% kMFF_RT_GeometryEGIG
% kMFF_RT_AnatomyEGIA

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function [EEG, data, events, chanlocs, mff] = mff_import(mffFile)

%matVer = ver('MATLAB');
%if datenum(matVer.Date) < 735595 % Matlab 2014a
%    error('This version of Matlab is too old. Use version 2014a or later');
%end

if nargin < 1 
    help mff_import;
    return;
end

% add full path if possible
mffPath = fileparts(mffFile);
if isempty(mffPath)
    mffFile = fullfile(pwd, mffFile);
elseif exist(fullfile(pwd, mffFile))
    mffFile = fullfile(pwd, mffFile);
end

% import data
[floatData, allDataSize, srate, nchans] = mff_importsignal(mffFile);

% aggregate all blocks
npns      = nchans(2);
nChannels = size(floatData{1}, 1);
blockSamples = cellfun(@(x)size(x,2), floatData);
if ~isequal(allDataSize/nChannels/4, blockSamples) && ...
        ~isequal(allDataSize/(nChannels+1)/4, blockSamples) % sometimes an empty PNS channel is removed
    error('Block sample size problem');
end
blockSamples = [0 cumsum(blockSamples)];
floatData = [ floatData{:} ];

if exist('eeg_emptyset.m', 'file')
    EEG = eeg_emptyset;
else
    EEG = [];
end
[~, EEG.setname] = fileparts(mffFile);
EEG.nbchan = double(nChannels);
EEG.data = floatData;
EEG.trials = 1;
EEG.srate = double(srate(1));
EEG.pnts  = size(EEG.data,2);
EEG.xmin  = 0;
EEG.xmax  = 1;
if exist('eeg_checkset.m', 'file') && exist('eeglab_options.m', 'file')
    EEG = eeg_checkset(EEG);
end

% scale signal with calibration values if necessary
disp('Make sure you to use the MFF ressource path name not the relative path (crash often related to that issue) ******');
info1 = mff_importinfon(mffFile, 1);
info2 = mff_importinfon(mffFile, 2);
if isfield(info1, 'calibration')
    disp('Calibrating data...');
    for iChan = 1:length(info1.calibration)
        floatData(iChan,:,:) = floatData(iChan,:,:)*info1.calibration(iChan);
    end
end
if isfield(info2, 'calibration')
    disp('Calibrating data...');
    for iChan = 1:length(info2.calibration)
        floatData(end-length(info2.calibration)+iChan,:,:) = floatData(end-length(info2.calibration)+iChan,:,:)*info2.calibration(iChan);
    end
end
info1.calibration = [];
info2.calibration = [];
EEG.etc.info1 = info1;
if ~isempty(info2) && length(fieldnames(info2)) > 1 
    EEG.etc.info2 = info2; % more than just calibration
end

% import info file
info    = mff_importinfo(mffFile);
layout  = mff_importsensorlayout(mffFile);
subject = mff_importsubject(mffFile);
begtime            = info.recordtimematlab;
EEG.etc.timezone   = info.timezone;
EEG.etc.mffversion = info.version;
EEG.etc.layout     = layout;
EEG.etc.subject    = subject;

% import coordinate layout
[EEG.chanlocs, EEG.ref] = mff_importcoordinates(mffFile);
if iscell(EEG.ref)
    EEG.ref = sprintf('%s ', EEG.ref{:});
end
EEG.urchanlocs = EEG.chanlocs;
pnschans                = mff_importpnsset(mffFile);
if length(pnschans) ~= npns && ~(length(pnschans) == size(EEG.data,1) && isempty(EEG.chanlocs))
    if length(pnschans) == npns+1
        % last PNS status channel missing because blank and removed by mff_importsignal
        pnschans(end) = [];
    else
        error('Number of PNS raw data channels is not equal to number of PNS channels');
    end
end
if ~isempty(EEG.chanlocs)
    if exist('eeg_checkchanlocs.m', 'file')
        EEG = eeg_checkchanlocs(EEG); % put fiducials in chanfinfo
        nChannels = length(EEG.chanlocs);
    end
end
if ~isempty(pnschans)
    if isempty(EEG.chanlocs) && length(pnschans) == size(EEG.data,1)
        % Only PNS channels
        EEG.chanlocs = pnschans;
    else
        if isempty(EEG.chanlocs)
            for iChan = 1:nChannels
                EEG.chanlocs(end+1).labels = [ 'E' num2str(iChan) ];
                EEG.chanlocs(end  ).type   = 'EEG';
            end
        end
        for iChan = 1:npns
            EEG.chanlocs(end+1).labels = pnschans(iChan).labels;
            EEG.chanlocs(end  ).type   = pnschans(iChan).type;
        end
    end
end
if exist('pop_chanedit', 'file')
    EEG=pop_chanedit(EEG, 'forcelocs',[],'nosedir','+Y');
else
    EEG.chaninfo.nosedir = '+Y';
end

EEG.etc.recordingtime = begtime;

% import events
[EEG.event, newtimezone] = mff_importevents(mffFile, begtime, EEG.srate);
if ~isequal(EEG.etc.timezone, newtimezone) && ~isempty(newtimezone)
    error('Time zone issue');
end
%mff_exportevents(EEG.event, 'test', EEG.etc.recordingtime, EEG.etc.timezone, EEG.srate);

% import continuous or epoch data
cont = mff_importepochs(mffFile, info.version);

% import continuous or epoch data
cat = mff_importcategories(mffFile, info.version);

% calculate epoch length
allEpochLen = [ [cont.endtime] - [cont.begintime] ];
if length(unique(allEpochLen)) > 1
    fprintf([ 'IMPORTANT Warning: cannot import trials of different length\n' ... 
              '  importing as segmented data (trial/category info will be lost)\n' ] );
    cat = [];
end

if ~isempty(cat)
    if sum(cellfun(@length, { cat.trials })) ~= length(cont)
        error('The number of "segments" do not match the number of epochs');
    end
    
    % check time consistency
    iTrial = ones(1, length(cat));
    epochLen = cont(1).endtime-cont(1).begintime;
    EEG.xmin = -(cat(1).trials(1).eventbegin-cat(1).trials(1).begintime)/1000000;
    epochDiffLat = cell(1, length(cat));
    EEG.pnts = round(epochLen/1000000*EEG.srate);
    lastOriEventIndex = length(EEG.event);
    
    % recorting categories
    catCont = [];
    for iCat = 1:length(cat)
        catContTmp = cat(iCat).trials;
        [catContTmp.name] = deal(cat(iCat).name);
        try
            catCont = [ catCont catContTmp ];
        catch
            fieldsTmp = fieldnames(catContTmp);
            for iCount = 1:length(catContTmp)
                for iField = 1:length(fieldsTmp)
                    catCont(end+1).(fieldsTmp{iField}) = catContTmp(iCount).(fieldsTmp{iField});
                end
            end
        end
    end
    [tmp, indices] = sort([catCont.begintime]);
    catCont = catCont(indices);
    
    for iBound = 1:length(cont) % do not add first event
        
        if catCont(iBound).begintime ~= cont(iBound).begintime
            fprintf('Warning: categornies and epoch information does not match for epoch\n', iBound);
        end
                
        % there is a natural jitter of a few millisecond for each time-locking
        % event within the uncertainty of the EEG sampling rate
        % epochDiffLat{iCat}(end+1) = cat(iCat).trials(iTrial(iCat)-1).eventbegin-cat(iCat).trials(iTrial(iCat)-1).begintime;
        if catCont(iBound).eventbegin-catCont(iBound).begintime ~= EEG.xmin
            % disp('Time locking event offset');
        end
        
        % check latency and block consistency
        sampleCalculated = (cont(iBound).begintime/1000000)*EEG.srate;
        sampleBlock      = blockSamples(cont(iBound).firstblock); % this assumes block of size 1
        if abs(sampleCalculated-sampleBlock) > 1e-10
            fprintf('Warning: segment discontinuity (%d samples missing - pause in the recording or bug?)\n', iBound, sampleCalculated-sampleBlock);
        end
        
        % adding new fields to event structure
        trial = catCont(iBound);
        EEG.event(end+1).type   = trial.name;
        EEG.event(end).latency  = -EEG.xmin*EEG.srate+1+EEG.pnts*(iBound-1);
        EEG.event(end).duration = (trial.eventend-trial.eventbegin)/1000000*EEG.srate; % this is sometimes off by 1e-13

        EEG.event(end).epoch    = iBound;
        
        % copy other fields
        trialFields = setdiff(fields(trial), { 'name', 'begintime', 'endtime', 'eventbegin', 'eventend' });
        for iField = 1:length(trialFields)
            EEG.event(end).(trialFields{iField}) = trial.(trialFields{iField});
        end

        cont(iBound).samplebeg = cont(iBound).begintime/1000000*EEG.srate;
        cont(iBound).sampleend = cont(iBound).endtime/1000000*EEG.srate;
        cont(iBound).samplelen = cont(iBound).sampleend - cont(iBound).samplebeg;
        
        %EEG.event(end).latency  = (trial.eventbegin)/1000000*EEG.srate; % this is sometimes off by 1e-13
        %EEG.event(end).duration = (trial.eventend-trial.eventbegin)/1000000*EEG.srate; % this is sometimes off by 1e-13
    end
    
    % adjust event latencies in case of gaps between segments
    if any([cont(2:end).begintime] - [cont(1:end-1).endtime])
        for iEvent = 1:lastOriEventIndex
            % find closest sample
            inds = find( EEG.event(iEvent).latency > [cont.samplebeg]);
            if length(inds) > 1
                correction = (cont(inds(end)).samplebeg - sum([cont(1:inds(end)-1).samplelen]));
                if correction
                    EEG.event(iEvent).latency = EEG.event(iEvent).latency - correction;
                end
            end
        end
    end
    
    %% scan events and assign epoch
    for iEvent = 1:length(EEG.event)
        newepoch = floor((EEG.event(iEvent).latency+0.000001)/EEG.pnts)+1; % adding 1/1000000 of a sample is necessary to prevent the error below for some files
        if ~isempty(EEG.event(iEvent).epoch) && ~(newepoch == EEG.event(iEvent).epoch)
            %error(sprintf('Event %d, wong epoch index %d vs %d\n', iEvent, newepoch, EEG.event(iEvent).epoch));
            % This could be a bug in writing the file or a feature
            % to fix this, line 275, recompule the latency of the onset of the epoch
            % based on the block offset (since the latency of the onset of the epoch is wrong due either
            % to a bug or a paus in the recording) then use this information to recompute the event latency
            % (subtracting the wrong event begin and adding back the right one)
        end
        EEG.event(iEvent).epoch = newepoch;
    end
else
    % add boundary events
    discontinuities = 0;
    if ~isempty(cont) && cont(1).begintime ~= 0
        error('First discontinuity does not start at time 0');
    end
    for iBound = 2:length(cont) % do not add first event
        if cont(iBound-1).endtime ~= cont(iBound).begintime
            EEG.event(end+1).type = 'boundary';
            eventDuration   = cont(iBound).begintime - cont(iBound-1).endtime; % in microseconds
            discontinuities = discontinuities + eventDuration; 
            EEG.event(end).duration = eventDuration/1000000*EEG.srate; % in samples
            
            sampleCalculated = ((cont(iBound).begintime-discontinuities)/1000000)*EEG.srate;
            sampleBlock      = blockSamples(cont(iBound).firstblock); % this assumes block of size 1
            if abs(sampleCalculated-sampleBlock) > 1e-10
                fprintf('Warning: segment discontinuity (%d samples missing - pause in the recording or bug?)\n', sampleCalculated-sampleBlock);
            end
            EEG.event(end).latency  = (cont(iBound).begintime)/1000000*EEG.srate; % absolute time allow resorting events later
%            EEG.event(end).latency  = sampleCalculated;
        end
    end
    
    % rename Break cnt events
end

%% resort events and check event structure
if ~isempty(EEG.event)
    if any(cellfun(@isempty, { EEG.event.latency }))
        error('Some empty event latency')
    end
    [~,iEvent] = sort([EEG.event.latency]);
    EEG.event = EEG.event(iEvent);
    
    % remove duration of remove data portions from events
    subLatency = 0;
    for iEvent = 1:length(EEG.event)
        if strcmpi(EEG.event(iEvent).type, 'boundary')
            subLatency = subLatency + EEG.event(iEvent).duration;
        end
        EEG.event(iEvent).latency = EEG.event(iEvent).latency-subLatency;
    end
    
    % events are not checked in this function if called from pop_mffimport.m
    % because they would alter time locking event values when importing
    % data trials
    s = dbstack;
    if length(s) <= 1 || ~strcmpi(s(2).file, 'pop_mffimport.m')
        if exist('eeg_checkset.m', 'file') && exist('eeglab_options.m', 'file')
            EEG = eeg_checkset(EEG,'eventconsistency');
        end
    end
end

data = EEG.data;
events = EEG.event;
chanlocs = EEG.chanlocs;
mff = EEG.etc;
try
    EEG = eeg_checkset(EEG); 
catch
    EEG.pnts   = size(EEG.data,2);
    EEG.trials = size(EEG.data,3);
end