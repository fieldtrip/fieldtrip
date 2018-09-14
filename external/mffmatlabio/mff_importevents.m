% mff_importevents - import information from MFF 'eventsxxxxx.xml'
%                    files. If several files are detected, information
%                    from all files is imported.
%
% Usage:
%   [events, timezone] = mff_importevents(mffFile);
%
% Inputs:
%  mffFile - filename/foldername for the MFF file
%
% Output:
%  events   - EEGLAB event structure
%  timezone - time zone

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

function [events, timeZone] = mff_importevents(mffFile, begTime, srate)

events = [];
timeZone = [];
p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
import com.egi.services.mff.api.MFFFactory;
import com.egi.services.mff.api.MFFResourceType;
import com.egi.services.mff.api.LocalMFFFactoryDelegate;
import com.egi.services.mff.utility.ResourceUnmarshalException;
import com.egi.services.mff.api.Signal;
import com.egi.services.mff.api.SignalBlock;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
warning('on', 'MATLAB:Java:DuplicateClass');

% create a factory
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

% Create Signal object and read in event track file.
% Note that 3 is the value associated with the event track resource type.
eventFile = dir( fullfile(mffFile, 'Events_*.xml'));
if isempty(eventFile), return; end
timeZone = [];
eventCount = 1;
showWarning = true;

% use vararg2str to save events so they may be exported
eeglabExport = true;
if ~exist('vararg2str', 'file')
    eeglabExport = false;
end

for iEvent = 1:length(eventFile)
    eventtrackfilename = fullfile(mffFile, eventFile(iEvent).name);
    eventtracktype = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_EventTrack'));
    eventtrackObj = mfffactory.openResourceAtURI(eventtrackfilename, eventtracktype);
    
    if eventtrackObj.loadResource()
        
        eventlist = eventtrackObj.getEvents();
        nevents = eventlist.size();
        fprintf('Importing %d events from file %s...\n', nevents, eventFile(iEvent).name);
        
        % Read through each event track and inspect count of events
        % and beginTime for last event in each.
        if nevents>0

            multiplier = 24*60*60*srate;
            events(nevents).description = ''; % create last event
            
            for eventnum = 1:nevents
                
                eventObj = eventlist.get(eventnum-1);
                
                % Get keys for event and display key codes
                events(eventCount).begintime    = char(eventObj.getBeginTime());
                if events(eventCount).begintime(end-1) == '-'
                    events(eventCount).begintime = [ events(eventCount).begintime(1:end-1) '0' events(eventCount).begintime(end) ':00' ];
                elseif events(eventCount).begintime(end-2) == '-'
                    events(eventCount).begintime = [ events(eventCount).begintime(1:end-2) events(eventCount).begintime(end) ':00' ];
                end
                events(eventCount).classid      = char(eventObj.getClassID());
                events(eventCount).code         = char(eventObj.getCode());
                events(eventCount).description  = char(eventObj.getDescription());
                events(eventCount).duration     = eventObj.getDuration();
                events(eventCount).label        = char(eventObj.getLabel());
                events(eventCount).relativebegintime = eventObj.getRelativeBeginTime();
                events(eventCount).sourcedevice = char(eventObj.getSourceDevice());
                
                % compute latency in days with ms -> convert to samples
                % eventCount = 1; 
                events(eventCount).latency = (mff_decodetime(events(eventCount).begintime)-begTime)*multiplier;
                
                % dual time recoding
%                 tmp = mff_encodetime(events(eventCount).latency/multiplier+begTime, '08:00')
%                 fprintf('%s\n%s\n', events(eventCount).begintime, tmp);
                
                events(eventCount).type    = events(eventCount).code;
                
                % import keys
                keylist = eventObj.getKeys();
                events(eventCount).mffkeys = char(keylist);
                eventkeycount = keylist.size;
                keyVals = [];
                for q = 0:eventkeycount-1
                    theKey = keylist.get(q);
                    keyVals(q+1).code = char(theKey.getCode);
                    keyVals(q+1).data = char(theKey.getData);
                    keyVals(q+1).datatype = char(theKey.getDataType);
                    keyVals(q+1).description = char(theKey.getDescription);
                    cleancode = keyVals(q+1).code;
                    cleancode( cleancode < 48 ) = []; % invalid char
                    cleancode( cleancode > 57 & cleancode < 64 ) = []; % invalid char
                    try
                        events(eventCount).( [ 'mffkey_' cleancode ]) = keyVals(q+1).data;
                    catch
                        if showWarning
                            disp('Warning: issue when converting MFF event key ************');
                            showWarning = false;
                        end
                    end
                end
                if eeglabExport
                    events(eventCount).mffkeysbackup = vararg2str(keyVals); % for exporting
                end
                
                eventCount = eventCount+1;
                
            end
        end
    else
        fprintf('Could not load event file %s\n', eventtrackfilename);
    end
end
if ~isempty(events)
    % get the time zone (duplicate code in mff_importevents and mff_importinfo)
    minusSign = find(events(end).begintime == '+');
    if isempty(minusSign)
        minusSign = find(events(end).begintime == '-');
        minusSign = minusSign(end);
    end
    timeZone = events(end).begintime(minusSign(end):end);
    if length(timeZone) > 6
        timeZone =  [];
        disp('Issue with decoding the time zone');
    end
end

% % same as above but 10% faster. Was not retained because the increased
% % complexitity. Does not support importing keys
% 
% begintime         = cell(1,eventcount);
% classid           = cell(1,eventcount);
% code              = cell(1,eventcount);
% description       = cell(1,eventcount);
% duration          = cell(1,eventcount);
% label             = cell(1,eventcount);
% relativebegintime = cell(1,eventcount);
% sourcedevice      = cell(1,eventcount);
% eventkeys         = cell(1,eventcount);
% for eventnum = 1:eventcount
% 
%     eventObj = eventlist.get(eventnum-1);
%     
%     % Get keys for event and display key codes
%     begintime{eventnum}    = eventObj.getBeginTime();
%     classid{eventnum}      = eventObj.getClassID();
%     code{eventnum}         = eventObj.getCode();
%     description{eventnum}  = eventObj.getDescription();
%     duration{eventnum}     = eventObj.getDuration();
%     label{eventnum}        = eventObj.getLabel();
%     relativebegintime{eventnum} = eventObj.getRelativeBeginTime();
%     sourcedevice{eventnum} = eventObj.getSourceDevice();
%     eventkeys{eventnum}    = eventObj.getKeys();
%     
% end
% 
% disp('Converting events to Matlab format...')
% begintime    = cellfun(@(x)char(x), begintime, 'uniformoutput', false);
% begintime2   = cellfun(@(x)mff_decodetime(x), begintime, 'uniformoutput', false);
% classid      = cellfun(@(x)char(x), classid, 'uniformoutput', false);
% code         = cellfun(@(x)char(x), code, 'uniformoutput', false);
% description  = cellfun(@(x)char(x), description, 'uniformoutput', false);
% label        = cellfun(@(x)char(x), label, 'uniformoutput', false);
% relativebegintime = cellfun(@(x)char(x), relativebegintime, 'uniformoutput', false);
% sourcedevice = cellfun(@(x)char(x), sourcedevice, 'uniformoutput', false);
% eventkeys    = cellfun(@(x)char(x), eventkeys, 'uniformoutput', false);
% 
% events = struct('begintime', begintime, 'latency', begintime2, 'classid', classid, 'code', code, 'description', description, ...
%                 'label', label, 'relativebegintime', relativebegintime, 'sourcedevice', sourcedevice, 'eventkeys', eventkeys);

% Save as above but not faster
% eventlist = eventlist.toArray;
% eventlist = cell(eventlist);
% 
% begintime    = cellfun(@(x)char(x.getBeginTime()), eventlist, 'uniformoutput', false);
% begintime2   = cellfun(@(x)mff_decodetime(x), begintime, 'uniformoutput', false);
% classid      = cellfun(@(x)char(x.getClassID()), eventlist, 'uniformoutput', false);
% code         = cellfun(@(x)char(x.getCode()), eventlist, 'uniformoutput', false);
% description  = cellfun(@(x)char(x.getDescription()), eventlist, 'uniformoutput', false);
% label        = cellfun(@(x)char(x.getLabel()), eventlist, 'uniformoutput', false);
% sourcedevice = cellfun(@(x)char(x.getSourceDevice()), eventlist, 'uniformoutput', false);
% eventkeys    = cellfun(@(x)char(x.getKeys()), eventlist, 'uniformoutput', false);
% 
% relativebegintime = cellfun(@(x)x.getRelativeBeginTime(), eventlist, 'uniformoutput', false);
% duration          = cellfun(@(x)x.getDuration(), eventlist, 'uniformoutput', false);
% 
% events = struct('begintime', begintime, 'latency', begintime2, 'classid', classid, 'code', code, 'type', code, 'description', description, ...
%                 'label', label, 'duration', duration, 'relativebegintime', relativebegintime, 'sourcedevice', sourcedevice, 'eventkeys', eventkeys);

