%% write_mff_event.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012
%  Copyright 2012 EGI. All rights reserved.
%  Writes an event file to the MFF.
%
%  This function has several limitations: 
%  -- It doesn't necessarily fit into the Field Trip framework, although it
%  the newData and hdr parameters are Field Trip style variables.  
%  -- It doesn't modify the MFF file's history data. 
%
%%
function write_mff_event(filePath, eventfilename, trackType, trackName, events, hdr)
if ~strcmp(trackType, 'EVNT') && ~strcmp(trackType, 'STIM') && ~strcmp(trackType, 'PAT ')
    fprintf('ERROR - third parameter (trackType) must be on of "EVNT", "STIM" or "PAT ".\n');
    return;
end
if isempty(hdr)
    summaryInfo = mff_getSummaryInfo(filePath);
else
    summaryInfo = hdr.orig;
end
infoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info, 'info.xml', filePath);
beginTimeStr = infoObj.getRecordTime();
MFFUtil = javaObject('com.egi.services.mff.utility.MFFUtil');
dstURI = [filePath filesep 'Events_' eventfilename '.xml'];

fields = fieldnames(events);
hasDurationRemainders = true;
if isempty(find(strcmp(fields, 'durationRemainder'), 1));
    hasDurationRemainders = false;
end
hasSampleRemainders = true;
if isempty(find(strcmp(fields, 'sampleRemainder'), 1));
    hasSampleRemainders = false;
end

%%
newEventList = javaObject('java.util.ArrayList')
numEvents = size(events,2)
for p = 1:numEvents
    if isempty(events(p).value)
        event = javaObject('com.egi.services.mff.api.Event');
        event.setCode(events(p).type);
        
        doDurationRemainder = false;
        if hasDurationRemainders
            if events(p).durationRemainder ~= 0
                doDurationRemainder = true;
            end
        end
        duration = events(p).duration;
        if doDurationRemainder
            duration = samples2Micros(duration - 1, summaryInfo.sampRate) + events(p).durationRemainder;
        else
            duration = samples2Micros(events(p).duration, summaryInfo.sampRate);
        end
        event.setDuration(duration);

        % get begin time in microsecs
        % convert from epoch sample (as if no time went by during breaks)
        % to sample of continuous. 
        epochSample = events(p).sample;
        sample = epochSample2Sample(epochSample, summaryInfo.epochBeginSamps, summaryInfo.epochNumSamps);
        % Convert from samples to microsecs
        microsecs = samples2Micros(sample, summaryInfo.sampRate);
        if hasSampleRemainders
            microsecs = microsecs + events(p).sampleRemainder;
        end

        % String timeInStringForm;
        % long timeMicrosecondForm;
        % long newTime = timeMicrosecondForm  + MFFUtil.getTimeInMicroseconds(timeInStringForm);
        beginTimeMicrosecs = uint64(MFFUtil.getTimeInMicroseconds(beginTimeStr));
        microsecsDate = uint64(microsecs + beginTimeMicrosecs);
        % String newTimeInStringForm = MFFUtil.getDateTime(newTime, null);
        % Note that the null argument will give you the default time zone. If you want a timezone from the source, then do:
        % String newTimeInStringForm = MFFUtil.getDateTime(newTime, MFFUtil.getTimeZone(timeInStringForm));
        dateTimeStr = MFFUtil.getDateTime(microsecsDate, MFFUtil.getTimeZone(beginTimeStr));
        event.setBeginTime(dateTimeStr);

        newEventList.add(event);
    end
end

% add events to event object

delegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
factory = javaObject('com.egi.services.mff.api.MFFFactory', delegate);
resourceVal = com.egi.services.mff.api.MFFResourceType.kMFF_RT_EventTrack;
resourceType = javaObject('com.egi.services.mff.api.MFFResourceType', resourceVal);
%    fprintf('%s %s\n', char(URI), char(resourceType));
factory.createResourceAtURI(dstURI, resourceType);
newEventTrackObj = factory.openResourceAtURI(dstURI, resourceType);
if ~isempty(newEventTrackObj)
    
    newEventTrackObj.setEvents(newEventList);
    newEventTrackObj.setTrackType(trackType);
    newEventTrackObj.setName(trackName);
    newEventTrackObj.saveResource();

else
    fprintf('Could not create event track.\n');
end

%newEventTrackObj = javaObject('com.egi.services.mff.api.EventTrack', true);
%newEventTrackObj.setEvents(newEventList);
%newEventTrackObj.setTrackType(trackType);
%newEventTrackObj.setName(trackName);
%newEventTrackObj = newEventTrackObj.marshal(newEventTrackObj, dstURI);

function micros = samples2Micros(samples, sampRate)
sampDuration = 1000000/sampRate;
micros = uint64(samples*sampDuration);

function sampleNum = epochSample2Sample(epochSampleNum, epochBeginSamps, epochNumSamps)
epoch = 1;
while (epochSampleNum > sum(epochNumSamps(1:epoch)))
    epoch = epoch + 1;
end
sampleNum = ((epochSampleNum - sum(epochNumSamps(1:epoch-1))) + epochBeginSamps(epoch)) - 1;
