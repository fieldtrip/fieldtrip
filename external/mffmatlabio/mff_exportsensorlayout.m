% mff_exportsensorlayout - export information from MFF 'sensorLayout.xml' file
%
% Usage:
%   mff_exportsensorlayout(EEG, mffFile);
%
% Inputs:
%  EEG     - EEGLAB structure
%  mffFile - filename/foldername for the MFF file

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

function mff_exportsensorlayout(EEG, mffFile)

if ~isfield(EEG.etc, 'layout') || isempty(EEG.etc.layout)
    return;
end

matStruct = EEG.etc.layout;

% check for removed channels and update structure
if length(matStruct.sensors) ~= length(EEG.chanlocs)
    % different number of channels - scan for removed channels
    chanlocs = EEG.chanlocs;
    if isempty(chanlocs) || ~isfield(chanlocs, 'type') || isempty(chanlocs(1).type)
        currenChans = 1:size(EEG.data,1);
    else
        % find EEG channels and convert them
        eegChans = cellfun(@(x)isequal(lower(x), 'eeg'), {chanlocs.type});
        try
            currenChans = cellfun(@(x)str2num(x(2:end)), {chanlocs(eegChans).labels});
        catch
            % this can fail if EGI channels are not E1, E2, etc...
            currenChans = 1:size(EEG.data,1);
        end    
    end
    
    for iSens = 1:length(matStruct.sensors)
        if ismember(matStruct.sensors(iSens).number, currenChans)
            matStruct.sensors(iSens).type = 0;
        elseif matStruct.sensors(iSens).type == 0
            matStruct.sensors(iSens).type = 2;
        end
    end
    
end

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
mfffactory.createResourceAtURI(sLayoutURI, sensorLayoutType);
layoutObj = mfffactory.openResourceAtURI(sLayoutURI, sensorLayoutType);

variables = ...
    { ...
    'Name'            'char'  {};
    'OriginalLayout'  'char'  {} ;
    'Sensors'         'array' { 'Name' 'char' {}; 'Number' 'real' {}; 'X' 'real' {}; 'Y' 'real' {}; 'Z' 'real' {}; 'Type' 'real' {}; 'Identifier' 'real' {} };
    'Threads'         'array' { 'First' 'real' {}; 'Second' 'real' {} };
    'TilingSets'      'array' { '' 'array' {} };
    'Neighbors'       'array' { 'ChannelNumber' 'real' {}; 'Neighbors' 'array' { } } };

layout = [];
layoutObj = mff_setobj(layoutObj, matStruct, variables);
layoutObj.saveResource();

