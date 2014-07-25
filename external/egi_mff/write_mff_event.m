%% write_mff_event.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012, 4/15/2014
%  Copyright 2012, 2014 EGI. All rights reserved.
% 
%  Writes an events structure (as described at
%  http://fieldtrip.fcdonders.nl/reference/ft_read_event) to an event-track
%  file in an existing MFF file (filePath).
% 
%  filePath ? The path to the .mff file. 
% 
%  trackName ? the name of the track as it will appear in Net Station.
% 
%  replace ? a boolean that indicates what the code should do if the track
%  already exists. If set to ?true?, the existing track will get
%  overwritten. If set to ?false?, the code will end with an error message.
% 
%  hdr ? FieldTrip header. You have the option of passing in the header, or
%  []. If you pass in the header, it pulls data out of it, rather than
%  recomputing them.
% 
%  This function doesn't modify the MFF file's history data. To do that,
%  call mff_write_history, described below.
%%
function write_mff_event(filePath, trackName, events, replace, hdr)
if isempty(hdr)
    summaryInfo = mff_getSummaryInfo(filePath);
else
    summaryInfo = hdr.orig;
end
infoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info, 'info.xml', filePath);
beginTimeStr = infoObj.getRecordTime();
MFFUtil = javaObject('com.egi.services.mff.utility.MFFUtil');
dstURI = [filePath filesep 'Events_' pathSafe(trackName) '.xml'];
if ~replace
    if exist(dstURI, 'file') == 2
        theException = MException('EGI_MFF:EVENTTRACK_EXISTS', 'Specified event track exists. If you want to overwrite, set 4th parameter to true.');
        throw(theException);
    end
end

hasDurationRemainders = false;
hasSampleRemainders = false;
fields = fieldnames(events(1));
if ~isempty(find(strcmp(fields, 'orig'), 1));
    if ~isempty(find(strcmp(fields, 'durationRemainder'), 1));
        hasDurationRemainders = true;
    end
    if ~isempty(find(strcmp(fields, 'sampleRemainder'), 1));
        hasSampleRemainders = true;
    end
else
end

%%
newEventList = javaObject('java.util.ArrayList');
numEvents = size(events,2);
for p = 1:numEvents
    if isempty(events(p).value)
        event = javaObject('com.egi.services.mff.api.Event');
        type = events(p).type;
        typeLen = size(type, 2);
        if typeLen > 4
            typeOrig = type;
            type = type(1:4);
            fprintf('*** Event type field ''%s'' is over 4 characters. Truncating to ''%s''. ***\n', typeOrig, type);
        elseif typeLen < 4
            typeOrig = type;
            for q = typeLen+1:4
                type = [type '_'];
            end
            fprintf('*** Event type field ''%s'' is under 4 characters. Padding to ''%s''. ***\n', typeOrig, type);
        end
        event.setCode(type);
        
        doDurationRemainder = false;
        if hasDurationRemainders
            if events(p).orig.durationRemainder ~= 0
                doDurationRemainder = true;
            end
        end
        duration = events(p).duration;
        if doDurationRemainder
            duration = samples2Micros(duration - 1, summaryInfo.sampRate) + events(p).orig.durationRemainder;
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
            microsecs = microsecs + events(p).orig.sampleRemainder;
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
    newEventTrackObj.setTrackType('EVNT');
    newEventTrackObj.setName(trackName);
    newEventTrackObj.saveResource();

else
    fprintf('Could not create event track.\n');
end

function micros = samples2Micros(samples, sampRate)
sampDuration = 1000000/sampRate;
micros = uint64(samples*sampDuration);

function sampleNum = epochSample2Sample(epochSampleNum, epochBeginSamps, epochNumSamps)
epoch = 1;
while (epochSampleNum > sum(epochNumSamps(1:epoch)))
    epoch = epoch + 1;
end
sampleNum = ((epochSampleNum - sum(epochNumSamps(1:epoch-1))) + epochBeginSamps(epoch)) - 1;

% Unsafe Chars: <space>:/\?%*|\"<>
function pathSafeFile = pathSafe(inFile)
unsafeInds = find(...
    inFile == ' ' |...
    inFile == ':' |...
    inFile == '/' |...
    inFile == '\' |...
    inFile == '?' |...
    inFile == '%' |...
    inFile == '*' |...
    inFile == '|' |...
    inFile == '\' |...
    inFile == '"' |...
    inFile == '<' |...
    inFile == '>'...
);
pathSafeFile = inFile;
pathSafeFile(unsafeInds) = '_';
