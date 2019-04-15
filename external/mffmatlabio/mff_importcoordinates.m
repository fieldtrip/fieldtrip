% mff_importcoordinates - import information from MFF 'coordinates.xml' file
%
% Usage:
%   [chanlocs, ref] = mff_importcoordinates(mffFile);
%
% Inputs:
%  mffFile - filename/foldername for the MFF file
%
% Output:
%  chanlocs - EEGLAB channel location structure
%  ref      - reference

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

function [chanlocs, reference, rVal] = mff_importcoordinates(mffFile)

chanlocs = [];
rVal = true;

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
warning('on', 'MATLAB:Java:DuplicateClass');

% Create an MFFFactory object.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

% Create Signal object and read in event track file.
% Note that 3 is the value associated with the event track resource type.
sLayoutURI = fullfile(mffFile, 'coordinates.xml');
sensorLayoutType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Coordinates'));

cLayout = mfffactory.openResourceAtURI(sLayoutURI, sensorLayoutType);

listfid  = [];
ref      = [];
reference = {};
if ~isempty(cLayout)
    
    %try
    
    if cLayout.loadResource()
        
        %  The information about the Coordinate Layout Resource.
        fprintf('Coordinate acquisition Time: %s\n', char(cLayout.getAcquisitionTime()));
        fprintf('Coordinate acquisition Method: %s\n', char(cLayout.getAcquisitionMethod()));
        
        %  Check to see if we have a default subject.
        %  fprintf('Default Subject: %d\n', cLayout.getDefaultSubject());
        
        %  Get the sesnors.
        sLayout =  cLayout.getSensorLayout();
        
        if ~isempty(sLayout )
            
            %  Information about the sensor layout.
            
            fprintf('Coordinate Layout Name: %s\n', char(sLayout.getName()));
            
            %  If Montage
            if ~isempty(sLayout.getOriginalLayout())
                fprintf('Cordinate Original Layout Name: %s\n', char(sLayout.getOriginalLayout()) );
            end
            sensors = sLayout.getSensors();
            
            if ~isempty(sensors )
                
                for iSensor = 1:sensors.size
                    sensor = sensors.get(iSensor-1);
                    description = char(sensor.getName());
                    chanlocs(iSensor).labels = [ 'E' num2str(sensor.getNumber()) ];
                    chanlocs(iSensor).description = description;
                    chanlocs(iSensor).X = sensor.getX();
                    chanlocs(iSensor).Y = sensor.getY();
                    chanlocs(iSensor).Z = sensor.getZ();
                    chanlocs(iSensor).identifier = sensor.getIdentifier();
                    sensorType = sensor.getType();
                    switch sensorType
                        case 0, chanlocs(iSensor).type = 'EEG';
                        case 1, chanlocs(iSensor).type = 'EEG'; ref = [ref iSensor]; reference = { reference{:} chanlocs(iSensor).labels };
                        case 2, chanlocs(iSensor).type = 'FID'; listfid = [ listfid iSensor ];
                    end
                end
                
            end
        else
            fprintf('Error: Coordinate Layout is null\n');
            rVal = false;
        end
        
    else
        fprintf('Error: could not load the coordinates resource.\n');
        rVal = false;
    end
    
end

if ~isempty(ref)
    for iSensor = 1:length(chanlocs)
        if strcmpi(chanlocs(iSensor).type, 'EEG'), chanlocs(iSensor).ref = [ chanlocs(ref).labels ]; end
    end
end

if ~isempty(chanlocs)
    chanlocs = convertlocs(chanlocs, 'cart2all');
end


