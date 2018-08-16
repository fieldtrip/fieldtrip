% mff_importpnsset - import information from MFF 'pnsSet.xml' file
%
% Usage:
%   chaninfo = mff_importpnsset(mffFile);
%
% Inputs:
%  mffFile - filename/foldername for the MFF file
%
% Output:
%  chaninfo - EEGLAB channel location structure

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

function [chaninfo, rVal] = mff_importpnsset(mffFile)

chanlocs = [];
rVal = true;

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
warning('on', 'MATLAB:Java:DuplicateClass');

% Create an MFFFactory object.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

% Create Signal object and read in event track file.
% Note that 3 is the value associated with the event track resource type.
pnsURI = fullfile(mffFile, 'pnsSet.xml');
pnsType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_PNSSet'));
pns = mfffactory.openResourceAtURI(pnsURI, pnsType);

chaninfo = [];
if ~isempty(pns)
    
    if pns.loadResource()
        
        sensors = pns.getPNSSensors();
        
        if ~isempty(sensors )
            
            for iSensor = 1:sensors.size
                sensor = sensors.get(iSensor-1);
                chaninfo(iSensor).labels = char(sensor.getName);
                chaninfo(iSensor).type   = 'PNS';
            end
            
        end
    else
        fprintf('Error: pnsSet ressource is null\n');
        rVal = false;
    end
    
else
    
    rVal = false;
    
end
