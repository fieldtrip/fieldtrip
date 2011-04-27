%% Trigger search
% This first section is typical example code to search the trigger lines
% for actual triggers. Specifying the trigger lines is optional, but if
% you don't, then please check the results of sqdgettriggers() to verify
% that spurious trigger lines weren't found.

clear
sqdfile = 'mySqdFile.sqd'

disp('Finding trigger lines and triggers...')
[triggerLines, triggerCounts, triggersList, fs] = sqdgettriggers(sqdfile,160:190)

%% Create & Save Epoched Data
% This second section is typical example code to use the detected triggers to
% break the data up into epochs/trials. In addition to the triggersList created
% above, the function sqdmakeepochs() requires the pre-trigger and
% post-trigger duration (both as number of samples). Optional are the base
% name of the epoch Matlab .mat files to create, and a flag (in the form of 
% a text string) that specifies what should be stored in the epoch files:
% the full epoched data, the average of epoched data, or both.

nPretriggerSamples = 100;
nPosttriggerSamples = 1000;
outputFilenameRoot = strrep(sqdfile,'.sqd','.epoch.mat');
epochAvgFlag = 'epoch+average';

[epochFiles, epochInfo] = sqdmakeepochs(sqdfile, triggersList, nPretriggerSamples, nPosttriggerSamples, outputFilenameRoot, epochAvgFlag);
disp(epochInfo)

% Save relevant output to a .mat file so this section need not be called again
save(outputFilenameRoot,'epochFiles', 'epochInfo')

%% Read Epoched Data and Analyze Results: 3 Scenarios
% This section would typically be in a separate file from the two sections
% above. This is where the data is actually analyzed and so is likely to be
% run many times with many modifications. The sections above would most
% likely only be run once.

% All scenarios below rely on this first section.

clear

sqdfile = 'mySqdFile.sqd'
epochHeaderFile = strrep(sqdfile,'.sqd','.epoch.mat');

% Load epoch header information and put into conveniently named variables
eh = load(epochHeaderFile); %eh = epoch header
disp(eh.epochInfo)

nChannels = eh.epochInfo.nChannels
nStimulusTypes = eh.epochInfo.nStimulusTypes
fs = eh.epochInfo.fs % Hz
Ts = 1000/fs % ms

nDuration = eh.epochInfo.nDuration % samples
nPretrigger = eh.epochInfo.nPretrigger % samples
nPosttrigger = eh.epochInfo.nPosttrigger % samples

tDuration = eh.epochInfo.tDuration % samples
tPretrigger = eh.epochInfo.tPretrigger % samples
tPosttrigger = eh.epochInfo.tPosttrigger % samples

t = (-tPretrigger+Ts:Ts:tPosttrigger);

% Read epoched and averaged data

epochFiles = eh.epochFiles;

%% Scenario 1
%
% Choose to only look at average data, which can be put into a 3 dimensional 
% matrix (time x channel x stimulus).
% All epochs/trials are used (no artifact rejection) and we assume
% that the average across epochs was pre-calculated during the epoching step.

avgData = nan(nDuration,nChannels,nStimulusTypes);
for mStim = 1:nStimulusTypes
    loadedData = load(epochFiles{mStim});
    avgData(:,:,mStim) = loadedData.averagedata;
end

% Baseline correct on pre-trigger responses
prestimAvgData = avgData(1:nPretrigger,:,:);
baseline = mean(prestimAvgData,1); % or mean(avgData,1) if baseline w.r.t. entire duration
avgDataBaselineCorrected = avgData - repmat(baseline,[nDuration 1 1]);

% Butterfly plot for each stimulus type
figure(1)
for mStim = 1:nStimulusTypes
    subplot(nStimulusTypes,1,mStim);plot(t,avgDataBaselineCorrected(:,:,mStim))
    xlim([-nPretrigger nPosttrigger])
end

% Plot RMS over channels for each stimulus type
figure(2)
for mStim = 1:nStimulusTypes
    subplot(nStimulusTypes,1,mStim);plot(t,std(avgDataBaselineCorrected(:,:,mStim),1,2))
    xlim([-nPretrigger nPosttrigger])
end


%% Scenario 2
%
% Instead analyze raw epoched data
% Here, individual trials can be checked, e.g., for artifacts.

artifactThresh = 1200; % pT
avgData = nan(nDuration,nChannels,nStimulusTypes);
for mStim = 1:nStimulusTypes
    loadedData = load(epochFiles{mStim});
    goodEpochs = (max(max(abs(loadedData.data))) <= artifactThresh);
    avgData(:,:,mStim) = mean(loadedData.data(:,:,goodEpochs),3);
    clear loadedData
end

% Baseline correct on pre-trigger responses
prestimAvgData = avgData(1:nPretrigger,:,:);
baseline = mean(prestimAvgData,1); % or mean(avgData,1) if baseline w.r.t. entire duration
avgDataBaselineCorrected = avgData - repmat(baseline,[nDuration 1 1]);

% Butterfly plot for each stimulus type
figure(3)
for mStim = 1:nStimulusTypes
    subplot(nStimulusTypes,1,mStim);plot(t,avgDataBaselineCorrected(:,:,mStim))
    xlim([-nPretrigger nPosttrigger])
end

% Plot RMS over channels for each stimulus type
figure(4)
for mStim = 1:nStimulusTypes
    subplot(nStimulusTypes,1,mStim);plot(t,std(avgDataBaselineCorrected(:,:,mStim),1,2))
    xlim([-nPretrigger nPosttrigger])
end

%% Scenario 3
%

% This does the same analysis as in Scenario 2, but all epochs are
% retained, instead of being thrown away, after each stimulus is processed.
%
% For some datasets, this may be too memory intensive, in which case 
% Scenario 2 is better. Paradoxically, it can be faster to read the data
% off the disk each time it is needed than to keep it in memory.

dataAsCell = cell(nStimulusTypes,1);
for mStim = 1:nStimulusTypes
    dataAsCell{mStim} = load(epochFiles{mStim});
end

artifactThresh = 1200; % pT
cleanDataAsCell = dataAsCell;
for mStim = 1:nStimulusTypes
    badEpochs = (max(max(abs(dataAsCell{mStim}.data))) > artifactThresh);
    cleanDataAsCell{mStim}.data(:,:,badEpochs)=[];
end

avgData = nan(nDuration,nChannels,nStimulusTypes);
for mStim = 1:nStimulusTypes
    avgData(:,:,mStim) = mean(cleanDataAsCell{mStim}.data,3);
end

% Baseline correct on pre-trigger responses
prestimAvgData = avgData(1:nPretrigger,:,:);
baseline = mean(prestimAvgData,1); % or mean(avgData,1) if baseline w.r.t. entire duration
avgDataBaselineCorrected = avgData - repmat(baseline,[nDuration 1 1]);

% Butterfly plot for each stimulus type
figure(5)
for mStim = 1:nStimulusTypes
    subplot(nStimulusTypes,1,mStim);plot(t,avgDataBaselineCorrected(:,:,mStim))
    xlim([-nPretrigger nPosttrigger])
end

% Plot RMS over channels for each stimulus type
figure(6)
for mStim = 1:nStimulusTypes
    subplot(nStimulusTypes,1,mStim);plot(t,std(avgDataBaselineCorrected(:,:,mStim),1,2))
    xlim([-nPretrigger nPosttrigger])
end

% If the number of epochs/stimulus is constant, they can put into a 4
% dimensional matrix (time x channel x epoch x stimulus), instead of a
% cell, but this may not be necessary (or even possible, if artifactual 
% trials are excluded).
%
% nEpochs =  dataAsCell{1}.nEpochs;
% data = nan(nDuration,nChannels,nEpochs,nStimulusTypes);
% for mStim = 1:nStimulusTypes
%      data(:,:,:,mStim) = dataAsCell{mStim}.data;
% end

