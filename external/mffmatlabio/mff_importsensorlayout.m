% mff_importsensorlayout - import information from MFF 'sensorLayout.xml' file
%
% Usage:
%   layout = mff_importsensorlayout(mffFile);
%
% Inputs:
%  mffFile - filename/foldername for the MFF file
%
% Output:
%  layout - Matlab structure containing informations contained in the MFF file.

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

function [layout, rVal] = mff_importsensorlayout(mffFile)

layout = [];
rVal = true;

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
warning('on', 'MATLAB:Java:DuplicateClass');

% Create an MFFFactory object.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

% Create Signal object and read in event track file.
sLayoutURI = fullfile(mffFile, 'sensorLayout.xml');
sensorLayoutType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_SensorLayout'));

sLayout = mfffactory.openResourceAtURI(sLayoutURI, sensorLayoutType);
% 
% layout = [];
% if ~isempty(sLayout)
% 
%     % Load the Sensor Layout resource.
%     if sLayout.loadResource()
%         layout.name           = char(sLayout.getName());
% 
%         % Get information about this layout.
%         fprintf('Layout Name: ', char(sLayout.getName()));
% 
%         if ~isempty(sLayout.getOriginalLayout)
%             fprintf('Original Layout Name: %s\n', char(sLayout.getOriginalLayout()));
%             layout.originallayout = char(sLayout.getOriginalLayout());
%         end
% 
%         sensors = sLayout.getSensors();
%         neighbors = sLayout.getNeighbors();
%         threads = sLayout.getThreads(); %Not CPU threads.
%         tilingSets = sLayout.getTilingSets();
% 
%         % Accessing the sensors.
%         if ~isempty(sensors)
% 
%             for iSensor = 1:sensors.size
%                 sensor = sensors.get(iSensor-1);
% 
%                 layout.sensor(iSensor).name = char(sensor.getName());
%                 layout.sensor(iSensor).number = sensor.getNumber();
%                 layout.sensor(iSensor).X = sensor.getX();
%                 layout.sensor(iSensor).Y = sensor.getY();
%                 layout.sensor(iSensor).Z = sensor.getZ();
%                 layout.sensor(iSensor).type = sensor.getType();
%             end
%         end
% 
%         % Sensor Layout threads are the visible connections between sensors on the net.
%         if ~isempty(threads)
% 
%             for iThread = 1:threads.size
%                 thread = threads.get(iThread-1);
%                 layout.thread(iThread).first  = thread.getFirst();
%                 layout.thread(iThread).second = thread.getSecond();
%             end
% 
%         end
% 
%         % Tiling Sets.
%         if ~isempty(tilingSets)
% 
%             for iTilingSet = 1:tilingSets.size
%                 tilingSet = tilingSets.get(iTilingSet-1);
%                 fprintf('Tiling Set:----------\n');
%                 for iSet = 1:tilingSet.size
%                     layout.tilingSet(iSet).value = tilingSet.get(iSet-1);
%                 end
%             end
% 
%         end
% 
%         fprintf('\n');
% 
%         % Neighbors.
%         if ~isempty(neighbors)
% 
%             for iNeighbor = 1:neighbors.size
%                 neighbor = neighbors.get(iNeighbor-1);
%                 layout.neighbors(iNeighbor).channelnumber = neighbor.getChannelNumber();
%                 channelNumbers =  neighbor.getNeighbors();
%                 for iChan = 1:channelNumbers.size
%                     layout.neighbors(iNeighbor).value(iChan) = channelNumbers.get(iChan-1);
%                 end
% 
%             end
% 
%         end
% 
%     else
% 
%         fprintf('Error: Could not load the Sensor Layout\n');
%         rVal = false;
% 
%     end
% 
% else
% 
%     fprintf('Error: The Sensor Layout resource does not exist\n');
%     rVal = false;
% end
% return

variables = ...
    { ...
    'Name'            'char'  {};
    'OriginalLayout'  'char'  {} ;
    'Sensors'         'array' { 'Name' 'char' {}; 'Number' 'real' {}; 'X' 'real' {}; 'Y' 'real' {}; 'Z' 'real' {}; 'Type' 'real' {}; 'Identifier' 'real' {} };
    'Threads'         'array' { 'First' 'real' {}; 'Second' 'real' {} };
    'TilingSets'      'array' { '' 'array' {} };
    'Neighbors'       'array' { 'ChannelNumber' 'real' {}; 'Neighbors' 'array' {} } };

layout = [];
if ~isempty(sLayout)
    try
        if sLayout.loadResource()
            layout = mff_getobj(sLayout, variables);
        end
    catch
        disp('Failed to load Layout ressource');
    end
end
