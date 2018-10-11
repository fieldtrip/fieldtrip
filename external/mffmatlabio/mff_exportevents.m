% mff_exportevents - export MFF EEG event from EEGLAB structure to 
%                    'Events_exported_from_EEGLAB.xml' file
% Usage:
%   mff_exportevents(EEG, mffFile);
%
% Inputs:
%  EEG     - EEGLAB structure
%  mffFile - filename/foldername for the MFF file (MFF file/folder must
%            already exist)

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

function mff_exportevents(EEG, mffFile)

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
warning('on', 'MATLAB:Java:DuplicateClass');

events   = EEG.event;
if isfield(EEG.etc, 'recordingtime')
     begTime = EEG.etc.recordingtime;
else begTime = now;
end
if isfield(EEG.etc, 'timezone')
     timeZone = EEG.etc.timezone;
else timeZone = '00:00';
end
srate    = EEG.srate;

% create a factory
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

% Create Signal object and read in event track file.
fileName = 'Events_exported_from_EEGLAB.xml';
eventtrackfilename = fullfile(mffFile, fileName);
eventtracktype = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_EventTrack'));
if mfffactory.createResourceAtURI(eventtrackfilename, eventtracktype)
    disp('Success at creating the event file');
else
    disp('Could not create event file');
end

fprintf('Exporting %d events...\n', length(events));
eventtrackObj = mfffactory.openResourceAtURI(eventtrackfilename, eventtracktype);

jList = javaObject('java.util.ArrayList');
addLatency = 0;
multiplier = 24*60*60*srate;
keyWarning = true;
for iEvent = 1:length(events)
    if strcmpi(events(iEvent).type, 'boundary') % do not export boundary events
        addLatency = addLatency + events(iEvent).duration;
    else
        eventObj = javaObject('com.egi.services.mff.api.Event');

        % Get keys for event and display key codes
        if isfield(events(iEvent), 'code')
            eventObj.setCode(        events(iEvent).code );
        else
            eventObj.setCode(num2str(events(iEvent).type));
        end
        if isfield(events(iEvent), 'description'),  eventObj.setDescription( events(iEvent).description ); end
        if isfield(events(iEvent), 'duration'   ),  eventObj.setDuration(    events(iEvent).duration ); end
        if isfield(events(iEvent), 'label'),        eventObj.setLabel(       events(iEvent).label ); end
        if isfield(events(iEvent), 'sourcedevice'), eventObj.setSourceDevice(events(iEvent).sourcedevice ); end
        
        if isfield(events, 'mffkeysbackup') && ~isempty(events(iEvent).mffkeysbackup)
            jListKeys = javaObject('java.util.ArrayList');
            tmpKeys = eval(events(iEvent).mffkeysbackup);
            for iKey = 1:length(tmpKeys);
                keyObj = javaObject('com.egi.services.mff.api.Key');
                keyObj.setCode(tmpKeys(iKey).code);
                keyObj.setData(tmpKeys(iKey).data);
                keyObj.setDataType(tmpKeys(iKey).datatype);
                keyObj.setDescription(tmpKeys(iKey).description);
                jListKeys.add(keyObj);
            end
            eventObj.setKeys(jListKeys);
        elseif keyWarning
            disp('Warning: no MFF backup key structure found - MFF keys not exported');
            disp('         import keys when loading the data to be able to export them');
            keyWarning = false;
        end
        
        realLatency = (events(iEvent).latency + addLatency)/multiplier+begTime;
        tmp = mff_encodetime(realLatency, timeZone);
        if isfield(events(iEvent), 'begintime')
            if ~isequal(events(iEvent).begintime, tmp) && ~isempty(events(iEvent).begintime)
                fprintf('Note: exported event time %d differ from original one %s vs %s\n', iEvent, events(iEvent).begintime, tmp);
            end
        end
        eventObj.setBeginTime(mff_encodetime(realLatency, timeZone));

        jList.add(eventObj);
    end
end
eventtrackObj.setEvents(jList);
eventtrackObj.saveResource();
