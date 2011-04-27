function [epochFiles, epochInfo] = sqdmakeepochs(sqdfile, triggersList, nPretriggerSamples, nPosttriggerSamples, outputFilenameRoot, epochAvgFlag)

%SQDMAKEEPOCHS epochs MEG data, and/or averages the epochs.
%
%    [epochFiles, epochInfo] = SQDMAKEEPOCHS(sqdfile, triggersList, nPretriggerSamples, nPosttriggerSamples, outFileRoot, epochAvgFlag);
%
%    INPUT variables:
%    sqdfile: the file holding the MEG data
%    triggersList: a cell array (each cell corresponds to a
%       different trigger line) in which each cell is a list of the trigger
%       events. The trigger events are labeled by sample number not time.
%    nPretriggerSamples: the number of pretrigger samples to include in the epoch
%    nPosttriggerSamples: the number of posttrigger samples to include in the epoch
%    outputFilenameRoot: the output file name root used to generate the
%       epochFiles (see below). The default is created using sqdfile,
%       removing the final '.sqd', and appending '-epochData'
%    epochAvgFlag: a string equal to one of the follwing values:
%       'epoch', 'average', or 'epoch+average', which is used to determine
%       whether the output files contain the epochs, the average of the
%       epochs, or both. The default is 'epoch+average'.
%
%    OUTPUT variables:
%    epochFiles: the files containging the epoched (and/or averaged) data.
%       The names of the epochFiles are of the form 'epochData01.mat',
%       'epochData02.mat', etc., if outputFilenameRoot is "epochData'.
%    epochInfo: a Matlab structure holding the information relevant to the epochs:
%       epochInfo.fs: the sample frequency, in Hz.
%       epochInfo.tPretrigger: duration of the pretrigger epoch section (in ms)
%       epochInfo.tPosttrigger: duration of the pretrigger epoch section (in ms)
%       epochInfo.tDuration: total duration of the epoch (sum of pre- and post
%           trigger durations) (in ms)
%       epochInfo.nPretrigger: duration of the pretrigger epoch section (in samples)
%       epochInfo.nPosttrigger: duration of the pretrigger epoch section (in samples)
%       epochInfo.nDuration: total duration of the epoch (sum of pre- and post
%           trigger durations) (in samples)

%
%    Version 0.9 beta 3
%    12 Feb 2010 
%    
%    by Jonathan Z. Simon

nTriggerLines = length(triggersList);

if ~exist('outputFilenameRoot','var')
    outputFilenameRoot = [];
end
if isempty(outputFilenameRoot)
    [pathstr,fname,fext] = fileparts(sqdfile);
    outputFilenameRoot = fullfile(pathstr,[fname '-epochData']);
end

if ~exist('epochAvgFlag','var')
    epochAvgFlag = [];
end
if isempty(epochAvgFlag)
    epochAvgFlag = 'epoch+average';
end
switch lower(epochAvgFlag)
    case 'epoch'
        doEpoch = true;
        doAvg = false;
    case 'average'
        doEpoch = false;
        doAvg = true;
    case 'epoch+average'
        doEpoch = true;
        doAvg = true;
    otherwise
        error('epochAvgFlag must be one of the following strings: ''epoch'', ''average'', or ''epoch+average''.')
end

nChannels = 157;
nDuration = nPretriggerSamples + nPosttriggerSamples;

info = sqdread(sqdfile,'info');
fs = info.SampleRate; % Hz
Ts = 1000/fs; % ms

tPretriggerSamples = nPretriggerSamples*Ts;
tPosttriggerSamples = nPosttriggerSamples*Ts;
tDuration = nDuration*Ts;

epochInfo = struct(...
    'nPretrigger', nPretriggerSamples, 'nPosttrigger', nPosttriggerSamples, 'nDuration', nDuration,...
    'tPretrigger', tPretriggerSamples, 'tPosttrigger', tPosttriggerSamples, 'tDuration', tDuration,...
    'fs', fs, 'nChannels', nChannels, 'nStimulusTypes', nTriggerLines);

epochFiles = cell(nTriggerLines,1);
for mTrigType = 1:nTriggerLines
    epochFiles{mTrigType} = [outputFilenameRoot num2str(mTrigType,'%02d') '.mat'];
end

dataChannels = 1:nChannels;

for mTrigType = 1:nTriggerLines
    nTriggers = length(triggersList{mTrigType});
    epochdataAll = nan(nDuration,nChannels,nTriggers);
    for mTrigLoc = 1:nTriggers
        epochdata = sqdread(sqdfile,'Samples',triggersList{mTrigType}(mTrigLoc)+[-nPretriggerSamples nPosttriggerSamples-1]);
        epochdataAll(:,:,mTrigLoc) = epochdata(:,dataChannels);
    end
    if doAvg
        averagedata = mean(epochdataAll,3);
    else
        averagedata = [];
    end
    if doEpoch
        data = epochdataAll;
    else
        data = [];
    end
    disp(['Creating file: ' epochFiles{mTrigType}])
    nEpochs = nTriggers;
    save(epochFiles{mTrigType},'data','averagedata','nPretriggerSamples', 'nPosttriggerSamples','nDuration','tPretriggerSamples', 'tPosttriggerSamples','tDuration','nEpochs');
end
    disp('Done creating files.')

