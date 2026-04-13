clear all
clc
folder= fullfile(pwd,'Data', filesep);
files = dir([folder '*.bvrh']);
headerFiles = {files.name};

for i=1:length(headerFiles)
    [hdr, EEGs] = eeg_loadbvrf(folder, headerFiles{i});
    [EEGs, com] = pop_loadbvrf(folder, headerFiles{i}, 'participantId', [], 'sampleInterval', [], 'channelIndx', []);
end
