% mff_exportcoordinates - export MFF 'coordinates.xml' file from EEGLAB structure
%
% Usage:
%   mff_exportcoordinates(EEG, mffFile);
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

function mff_exportcoordinates(EEG, mffFile)

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
import com.egi.services.mff.api.MFFFactory;
import com.egi.services.mff.api.MFFResourceType;
import com.egi.services.mff.api.LocalMFFFactoryDelegate;
import com.egi.services.mff.utility.ResourceUnmarshalException;
import com.egi.services.mff.api.SensorLayout;
import com.egi.services.mff.api.Sensor;
import com.egi.services.mff.api.Key;
import com.egi.services.mff.api.Neighbor;
import java.util.ArrayList;
warning('on', 'MATLAB:Java:DuplicateClass');

% remove PNS channels
if ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs, 'type')
    allTypes = { EEG.chanlocs.type };
    allTypes = cellfun(@(x)num2str(x), allTypes, 'uniformoutput', false);
    pnsChans = strmatch('PNS', allTypes, 'exact')';
    EEG.chanlocs(pnsChans) = [];
end

if isempty(EEG.chanlocs) || ~isfield(EEG.chanlocs(1), 'X') || isempty(EEG.chanlocs(1).X)
    return;
end
EEG=pop_chanedit(EEG, 'forcelocs',[],'nosedir','-Y');
[~, ~, chanlocs]= eeg_checkchanlocs( EEG.chanlocs, EEG.chaninfo);

% create a factory
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

% Create Signal object and read in event track file.
% Note that 3 is the value associated with the event track resource type
fileName = 'coordinates.xml';
coordinatefilename = fullfile(mffFile, fileName);
coordinatetype = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Coordinates'));
if mfffactory.createResourceAtURI(coordinatefilename, coordinatetype)
    disp('Success at creating the coordinate file');
else
    disp('File already exist, cannot create coordinate file');
end

fprintf('Exporting %d coordinate...\n', length(chanlocs));
coordinateObj = mfffactory.openResourceAtURI(coordinatefilename, coordinatetype);

layoutObj = javaObject('com.egi.services.mff.api.SensorLayout');
layoutObj.setName('Exported from EEGLAB');

jList = java.util.ArrayList;
for iSensor = 1:length(chanlocs)
    sensorObj = javaObject('com.egi.services.mff.api.Sensor');
    
    sensorObj.setNumber(iSensor);
    if isfield(chanlocs, 'description') && ~isempty(chanlocs(iSensor).description)
        sensorObj.setName(chanlocs(iSensor).description);
    else
        sensorObj.setName(chanlocs(iSensor).labels);
    end
    sensorObj.setX(chanlocs(iSensor).X);
    sensorObj.setY(chanlocs(iSensor).Y);
    sensorObj.setZ(chanlocs(iSensor).Z);
    if isfield(chanlocs, 'identifier') && ~isempty(chanlocs(iSensor).identifier) && (chanlocs(iSensor).identifier ~= 0)
        sensorObj.setIdentifier(chanlocs(iSensor).identifier);
    end
    if ~isfield(chanlocs, 'type') || isempty(chanlocs(iSensor).type)
        chanlocs(iSensor).type = 0;
    end
    switch chanlocs(iSensor).type
        case 'EEG', sensorObj.setType(0);
        case 'FID', sensorObj.setType(2);
    end
    if strcmpi(chanlocs(iSensor).labels, EEG.ref)
        sensorObj.setType(1);
    end
    jList.add(sensorObj);
    
end
layoutObj.setSensors(jList);
coordinateObj.setSensorLayout(layoutObj);
coordinateObj.setAcquisitionMethod('Exported from EEGLAB');
try
    coordinateObj.setAcquisitionTime(mff_encodetime(EEG.etc.timezone, EEG.etc.recordingtime));
catch
    coordinateObj.setAcquisitionTime('2006-01-01T00:00:00.000000-08:00');
end
coordinateObj.saveResource();

