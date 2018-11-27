% mff_exportepochs - export MFF 'epochs.xml' file from EEGLAB structure
%
% Usage:
%   mff_exportepochs(EEG, mffFile);
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

function mff_exportepochs(EEG, mffFile)

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
warning('on', 'MATLAB:Java:DuplicateClass');

% Create a factory.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);

%% create Segment to load time
epochsRType = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Epochs'));
catURI = [ mffFile filesep 'epochs.xml' ];
jList = javaObject('java.util.ArrayList');
if mfffactory.createResourceAtURI(catURI, epochsRType)
    fprintf('Epochs.xml file created successfully\n');
else
    fprintf('Epochs.xml ressource already exist, overwriting\n');
end
epochResource = mfffactorydelegate.openResourceAtURI(catURI, epochsRType);

% continuous data: export each portion of data
% epoched data: eport each epoch

% get offsets
durations = [];
if EEG.trials == 1
    samples = [];
    if ~isempty(EEG.event) && isfield(EEG.event, 'type') && isstr(EEG.event(1).type)
        boundaryEvent = strmatch( 'boundary', { EEG.event.type }, 'exact');
        samples       = [ EEG.event(boundaryEvent).latency ];
    end
    samples = [ 0 samples EEG.pnts ];
    if isfield(EEG.event, 'duration')
        durations = [ EEG.event(boundaryEvent).duration ];
        durations = cumsum(durations);
        durations = [0 durations durations(end) ];
    else
        durations = zeros(1,length(samples));
    end
else
    samples   = EEG.pnts*[0:EEG.trials];
    durations = zeros(1,length(samples));
end

% write data
block = 1;
warnOnce = false;
epochLen = [];
for iSample = 2:length(samples)
    
    epochObj = javaObject('com.egi.services.mff.api.Epoch');
    
    % latency of block
    epochObj.setBeginTime( round((samples(iSample-1)+durations(iSample-1))/round(EEG.srate)*1000000) ); % in microsec
    epochObj.setEndTime( round((samples(iSample)+durations(iSample-1))/round(EEG.srate)*1000000) );
    
    % check epoch length
    epochLenTmp = epochObj.getBeginTime() - epochObj.getEndTime();
    if isempty(epochLen)
        epochLen = epochLenTmp;
    elseif epochLen ~= epochLenTmp
        if warnOnce
            disp('Warning: epoch length differ (normal for continuous data)');
            warnOnce = true;
        end
    end
 
    % length of block
    blockinc = ceil((samples(iSample)-samples(iSample-1))/65536);
    epochObj.setFirstBlock(block);
    epochObj.setLastBlock(block+blockinc-1);
    block = block+blockinc;
    
    jList.add(epochObj);
end

epochResource.setEpochs(jList); 
epochResource.saveResource();

