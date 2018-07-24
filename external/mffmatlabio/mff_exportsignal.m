% mff_exportsignal - export MFF EEG data from EEGLAB structure into 
%                    'signal1.bin'
%
% Usage:
%   mff_exportsignal(EEG, mffFile);
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

function mff_exportsignal(EEG, mffFile)

% find PNS channels
if ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs, 'type')
    allTypes = { EEG.chanlocs.type };
    allTypes = cellfun(@(x)num2str(x), allTypes, 'uniformoutput', false);
    pnsChans = strmatch('PNS', allTypes, 'exact')';
else
    pnsChans = [];
end

p = fileparts(which('mff_importsignal.m'));
warning('off', 'MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(p, 'MFF-1.2.2-jar-with-dependencies.jar'));
import com.egi.services.mff.api.*;
import com.egi.services.mff.utility.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.util.ArrayList;
warning('on', 'MATLAB:Java:DuplicateClass');

% Create a factory.
mfffactorydelegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
mfffactory = javaObject('com.egi.services.mff.api.MFFFactory', mfffactorydelegate);
binfilename = { [mffFile filesep 'signal1.bin'] [mffFile filesep 'signal2.bin'] };
if isempty(pnsChans) 
    binfilename(2) = [];
    chanRange{1} = [1:EEG.nbchan];
else
    chanRange{1} = setdiff([1:EEG.nbchan], pnsChans);
    chanRange{2} = pnsChans;
end 

% PNS signal only -> move to file signal1.bin
if isempty(chanRange{1})
    chanRange = { chanRange{2} [] };
    binfilename(2) = [];
end
    
% Create Signal object and read in signal1.bin file.
signalresourcetype = javaObject('com.egi.services.mff.api.MFFResourceType', javaMethod('valueOf', 'com.egi.services.mff.api.MFFResourceType$MFFResourceTypes', 'kMFF_RT_Signal'));

for iFile = 1:length(binfilename)
    if mfffactory.createResourceAtURI(binfilename{iFile}, signalresourcetype)
        fprintf('Signal binary file created successfully\n');
    else
        fprintf('Signal binary ressource already exist, overwriting\n');
    end

    signalResource = mfffactory.openResourceAtURI(binfilename{iFile}, signalresourcetype);
    jList = java.util.ArrayList;
    if round(EEG.srate) ~= EEG.srate
        fprintf('Warning: sampling frequency need to be rounded from %1.2f to %1.0f\n', EEG.srate, round(EEG.srate));
    end
    if ~isempty(signalResource)

        % continuous data: export each portion of data
        % epoched data: eport each epoch
        nChans = length(chanRange{iFile});
        
        % get offsets
        if EEG.trials == 1
            samples = [];
            if ~isempty(EEG.event) && isfield(EEG.event, 'type') && isstr(EEG.event(1).type)
                boundaryEvent = strmatch( 'boundary', { EEG.event.type }, 'exact');
                samples       = [ EEG.event(boundaryEvent).latency ];
            end
            samples = round([ 0 samples EEG.pnts ]); % in rare cases, not rounding generates numerical instabilities (file MMVTD_Continuous_EEG.mff)
        else
            samples = EEG.pnts*[0:EEG.trials];
        end

        % add additional blocks if the data block is too large
        newSamples = [];
        for iSample = 1:length(samples)-1
            len        = samples(iSample+1)-samples(iSample);
            tmpSamples = 0:65536:len;
            newSamples = [ newSamples samples(iSample)+tmpSamples ];
        end
        samples = [newSamples samples(end)];

        % write data
        for iSample = 2:length(samples)

            newBlock = javaObject('com.egi.services.mff.api.SignalBlock');

            newBlock.version          = 1;
            newBlock.dataBlockSize    = (samples(iSample)-samples(iSample-1))*nChans*4;
            %fprintf('Data block size: %d\n', newBlock.dataBlockSize);
            newBlock.numberOfSignals  = nChans;
            newBlock.headerSize       = 20+nChans*8; % this was calculated heuristically

            newBlock.offsets         = [0:nChans-1]'*(samples(iSample)-samples(iSample-1))*4; % int32
            newBlock.signalDepth     = int32(ones(nChans, 1)*32); % int32 or 4 bytes
            newBlock.signalFrequency = int32(ones(nChans, 1)*round(EEG.srate)); % int32

            data = single(EEG.data(chanRange{iFile}, samples(iSample-1)+1:samples(iSample)))';
            newBlock.data = typecast(data(:), 'int8'); % LITTLE ENDIAN HERE - BIG ENDIAN MIGHT BE A PROBLEM

            jList.add(newBlock);
        end
        signalResource.setNumberOfBlocks(length(samples)-1);
        signalResource.setSignalBlocks(jList);
        signalResource.saveResource();

    else
        error('Error: Can not open the signal resource.\n');
    end
end
