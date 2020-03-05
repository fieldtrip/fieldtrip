% mff_exportpnsset - import information from MFF 'pnsSet.xml' file
%
% Usage:
%   chaninfo = mff_importpnsset(EEG, mffFile)
%
% EEG     - EEGLAB structure
% mffFile - filename/foldername for the MFF file (MFF file/folder must
%           already exist)

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

function mff_exportpnsset(EEG, mffFile)

% find PNS channels
if isempty(EEG.chanlocs), return; end
if ~isfield(EEG.chanlocs, 'type'), return; end
allTypes = { EEG.chanlocs.type };
allTypes = cellfun(@(x)num2str(x), allTypes, 'uniformoutput', false);
pnsChans = strmatch('pns', lower(allTypes), 'exact')';
if isempty(pnsChans), return; end

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
warning('on', 'MATLAB:Java:DuplicateClass');

% Create an MFFFactory object.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

% Create Signal object and read in event track file.
% Note that 3 is the value associated with the event track resource type.
pnsURI = fullfile(mffFile, 'pnsSet.xml');
pnsType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_PNSSet'));
if mfffactory.createResourceAtURI(pnsURI, pnsType)
    disp('Success at creating the PNS sensor file');
else
    disp('File already exist, not create coordinate file');
end

pns = mfffactory.openResourceAtURI(pnsURI, pnsType);

if ~isempty(pnsChans)
            
    jList = javaObject('java.util.ArrayList');

    for iChan = 1:length(pnsChans)
        sensorObj = javaObject('com.egi.services.mff.api.PNSSensor');
        sensorObj.setName(EEG.chanlocs(pnsChans(iChan)).labels);
        sensorObj.setNumber(iChan);
        jList.add(sensorObj);
    end

    pns.setPNSSensors(jList);
    pns.saveResource();
    
else
    
    fprintf('Error: could not load the pnsSet.xml resource.\n');
    
end
