%% read_mff_event.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012, 4/15/2014
%  Copyright 2012, 2014 EGI. All rights reserved.
%
%  Takes the path to the data and returns the events in the structure
%  described at http://www.fieldtriptoolbox.org/reference/ft_read_event. 
%
%  filePath ? The path to the .mff file. 
%
%  hdr ? FieldTrip header. You have the option of passing in the header, or
%  []. If you pass in the header, it pulls data out of it, rather than
%  recomputing them.
%%
function events = read_mff_event(filePath, hdr)
if isempty(hdr)
    summaryInfo = mff_getSummaryInfo(filePath);
else
    summaryInfo = hdr.orig;
end
% Pull the information about the epochs out of the summary info
epochBeginSamps = summaryInfo.epochBeginSamps;
epochNumSamps = summaryInfo.epochNumSamps;
epochFirstBlocks = summaryInfo.epochFirstBlocks;
epochLastBlocks = summaryInfo.epochLastBlocks;

% Get the start time of the recording
infoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info, 'info.xml', filePath);
beginTime = infoObj.getRecordTime();

events = [];
eventInd = 0;
% Create the meta data events, ie epoch breaks and time 0 events. 
for p = 1:size(summaryInfo.epochBeginSamps,2)
    eventInd = eventInd + 1;
    events(eventInd).type = ['break ' summaryInfo.epochType];
    events(eventInd).sample = samples2EpochSample(summaryInfo.epochBeginSamps(p), epochBeginSamps, epochNumSamps);
    events(eventInd).value = summaryInfo.epochLabels{p};
    events(eventInd).offset = [];
    events(eventInd).duration = summaryInfo.epochNumSamps(p);
    events(eventInd).timestamp = [];
    events(eventInd).orig.sampleRemainder = 0;
    events(eventInd).orig.durationRemainder = 0;
    events(eventInd).orig.trackname = 'metadata';
    events(eventInd).orig.keys = {};
    eventInds(eventInd,1) = events(eventInd).sample;
    eventInds(eventInd,2) = eventInd;
    if strcmp(summaryInfo.epochType, 'var')
        eventInd = eventInd + 1;
        events(eventInd).type = 't0';
        events(eventInd).sample = samples2EpochSample((summaryInfo.epochTime0(p) + summaryInfo.epochBeginSamps(p)) - 1, epochBeginSamps, epochNumSamps);
        events(eventInd).value = 't0';
        events(eventInd).offset = [];
        events(eventInd).duration = 1;
        events(eventInd).timestamp = []; % or calculate string
        events(eventInd).orig.sampleRemainder = 0;
        events(eventInd).orig.durationRemainder = 0;
        events(eventInd).orig.trackname = 'metadata';
        events(eventInd).orig.keys = {};
        eventInds(eventInd,1) = events(eventInd).sample;
        eventInds(eventInd,2) = eventInd;
    end
end
% Go through all the event tracks, if any...
eventtracknamelist = summaryInfo.javaObjs.mfffileObj.getEventTrackList(false);
eventtrackcount = eventtracknamelist.size();
if eventtrackcount > 0
    % MFFUtil is used to for operations on string-based timestamp
    MFFUtil = javaObject('com.egi.services.mff.utility.MFFUtil');
    for tracknum = 0:eventtrackcount-1
        trackname = eventtracknamelist.get(tracknum);
        tracknameL = lower(trackname);
        if strcmp(tracknameL(end-3:end), '.xml')
            trackname = trackname(1:end-4);
            eventTrackObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_EventTrack, tracknameL, filePath);
            eventList = eventTrackObj.getEvents;
            numEvents = eventList.size;
            for p = 0:numEvents-1
                % Java arrays are 0 based
                theEvent = eventList.get(p);
                eventTime = theEvent.getBeginTime;

                %Need to convert to samples (get sampling rate up above)
                eventTimeInMicros = uint64(MFFUtil.getTimeDifferenceInMicroseconds(eventTime , beginTime));
                [eventTimeInSamples, sampleRemainder] = mff_micros2Sample(eventTimeInMicros, summaryInfo.sampRate);
                eventTimeInEpochSamples = samples2EpochSample(eventTimeInSamples, epochBeginSamps, epochNumSamps);
                if eventTimeInEpochSamples < 1
                    % Error: invalid sample number
                else
                    eventInd = eventInd + 1;
    %                fprintf('%d %d %d\n', eventTimeInSamples, eventTimeInEpochSamples, eventTimeInSamples-eventTimeInEpochSamples);
                    % Matlab arrays are 1 based
                    events(eventInd).type = char(theEvent.getCode);
                    events(eventInd).sample = eventTimeInEpochSamples;
                    events(eventInd).orig.sampleRemainder = sampleRemainder;
                    events(eventInd).value = [];
                    events(eventInd).offset = [];
                    [events(eventInd).duration, events(eventInd).orig.durationRemainder] = mff_micros2Sample(theEvent.getDuration, summaryInfo.sampRate);
                    if events(eventInd).orig.durationRemainder > 0
                        events(eventInd).duration = events(eventInd).duration + 1;
                    end
                    events(eventInd).timestamp = []; %eventTime;
                    events(eventInd).orig.trackname = trackname;
                    % New code for keys. Java arrays are 0 based. 
                    keylist = theEvent.getKeys;
                    eventkeycount = keylist.size;
%                     events(eventInd).keys = cell(eventkeycount, 4);
                    % Create an empty keylist in case there are no keys,
                    % e.g. DINs. 
                    events(eventInd).orig.keys = {};
                    for q = 0:eventkeycount-1
                        theKey = keylist.get(q);
                        theKeyCode = char(theKey.getCode);
                        theKeyData = char(theKey.getData);
                        theKeyDataType = char(theKey.getDataType);
                        theKeyDescription = char(theKey.getDescription);
                        events(eventInd).orig.keys{q+1}.code = theKeyCode;
                        events(eventInd).orig.keys{q+1}.data = theKeyData;
                        events(eventInd).orig.keys{q+1}.datatype = theKeyDataType;
                        events(eventInd).orig.keys{q+1}.description = theKeyDescription;
%                         events(eventInd).keys(q+1).code = theKeyCode;
%                         events(eventInd).keys(q+1).data = theKeyData;
%                         events(eventInd).keys(q+1).datatype = theKeyDataType;
%                         events(eventInd).keys(q+1).description = theKeyDescription;
%                         events(eventInd).keys{q+1, 1} = theKeyCode;
%                         events(eventInd).keys{q+1, 2} = theKeyData;
%                         events(eventInd).keys{q+1, 3} = theKeyDataType;
%                         events(eventInd).keys{q+1, 4} = theKeyDescription;
                    end
                    eventInds(eventInd,1) = events(eventInd).sample;
                    eventInds(eventInd,2) = eventInd;
                end
            end
        end
    end
    % Sort events, which includes metadata events and events that have come
    % from different tracks, by time. 
    eventInds = sortrows(eventInds);
    for p = 1:eventInd
        nextEventInd = eventInds(p,2);
        sortedEvents(p) = events(nextEventInd);
    end
    events = sortedEvents;
end

% Converts from samples since start of recording (as if there were no
% breaks) to samples in file. 
function epochSampleNum = samples2EpochSample(sampleNum, epochBeginSamps, epochNumSamps)
numEpochs = size(epochBeginSamps,2);
p = 1;
epochSampleNum = 0;
while (sampleNum > (epochBeginSamps(p) + epochNumSamps(p) - 1)) && (p < numEpochs)
    epochSampleNum = epochSampleNum + epochNumSamps(p);
    p = p+1;
end
if p <= numEpochs
    if sampleNum >= epochBeginSamps(p)
        epochSampleNum = epochSampleNum + ((sampleNum - epochBeginSamps(p))+1);
    else
        epochSampleNum = -1;
        % Error: sample falls between epochs
    end
else
    epochSampleNum = -2;
    % Error: sample is after last epoch
end



