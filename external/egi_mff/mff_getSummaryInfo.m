%% mff_getSummaryInfo.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012
%  Copyright 2012 EGI. All rights reserved.
%  Support routine for MFF Matlab code. Not intended to be called directly.
%  
%  Returns data about the MFF file that is used by the multiple functions
%  that read and write MFF files. Some of this data goes directly into the
%  field trip header. It pulls the data from the following MFF items.
% 
%  Each MFF has one or more bin files, where only one EEG. 
% 
%  Each EEG bin file has 1 or more blocks objects. 
% 
%  Blocks is an array of the block objects for a given bin file. 
% 
%  Each block object contains info about the sampling rate, number of
%  channels, and the actual voltage data. It is assumed that the sampling
%  rates and numbers of channels are the same for all blocks.
% 
%  Each MFF has an epochs file that has an array of epochs. 
% 
%  Each epoch has a begin and end time (in nanoseconds since recording
%  start), and first and last block.
% 
%  If the data have been segmented, then there is a categories file that has
%  an array of categories. Each category has an array of segments. There is
%  a 1-1 correspondence between segments and epochs. 
% 
%  Each segment has a begin time and end time in nanoseconds, which can be
%  used to map the segment to the corresponding epoch. Each segment also has
%  a time zero event with its time in microseconds since recording
%  start. 
%%
function summaryInfo = mff_getSummaryInfo(filePath)
% create the MFFFile object

mfffileObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_MFFFile, [], filePath);
%mfffileObj = javaObject('com.egi.services.mff.api.MFFFile', filePath, true);

[binObj blocks] = getEEGBlocks(mfffileObj, filePath);
numblocks = binObj.getNumberOfBlocks();
if numblocks > 0 % Should always be

    summaryInfo.mfffileObj = mfffileObj;
    summaryInfo.binObj = binObj; % bin obj for the EEG data
    summaryInfo.blocks = blocks; % block objects for the EEG bin obj

    blockObj = blocks.get(0); %zero based
    sampRate = double(blockObj.signalFrequency(1)); % 1 based
    nChans = blockObj.numberOfSignals;
    summaryInfo.sampRate = sampRate;
    summaryInfo.nChans = nChans;

    [epochType epochBeginSamps epochNumSamps epochFirstBlocks epochLastBlocks epochLabels epochTime0] = getEpochInfos(filePath, sampRate);
    summaryInfo.epochType = epochType;
    summaryInfo.epochBeginSamps = epochBeginSamps;
    summaryInfo.epochNumSamps = epochNumSamps;
    summaryInfo.epochFirstBlocks = epochFirstBlocks;
    summaryInfo.epochLastBlocks = epochLastBlocks;
    summaryInfo.epochLabels = epochLabels;
    summaryInfo.epochTime0 = epochTime0;

    % Check that number of channels and sampling rate is consistent
    % across blocks. 
    for x = 0:numblocks-1
        blockObj = blocks.get(x);
        sampRate = double(blockObj.signalFrequency(1)); % 1 based
        nChans = blockObj.numberOfSignals;
        if sampRate ~= summaryInfo.sampRate
            % Error: Inconsistent sampling rate. Should never occur.
            % todo?: error handling
        end
        if nChans ~= summaryInfo.nChans
            % Error: Inconsistent number of channels. Should never
            % occur. todo?: error handling
        end
    end
else
    % Error: Signal has 0 blocks. Should never occur. todo?: error
    % handling
end

% Returns the bin object corresponding to the EEG, and the blocks array
% associated associated with it. 
function [binObj blocks] = getEEGBlocks(mfffileObj, filePath)
EEGFile = mff_getEEGFilename(mfffileObj);
binObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal, EEGFile, filePath);
blocks = binObj.getSignalBlocks();

% Each return value is an array with one element for each epoch, except epochType, which is global description. 
function [epochType epochBeginSamps epochNumSamps epochFirstBlocks epochLastBlocks epochLabels epochTime0] = getEpochInfos(filePath, sampRate)
epochList = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Epochs, 'epochs.xml', filePath);
epochListArray = epochList.getEpochs;
numEpochs = epochListArray.size;
epochBeginSamps = zeros(1,numEpochs);
epochNumSamps = zeros(1,numEpochs);
epochFirstBlocks = zeros(1,numEpochs);
epochLastBlocks = zeros(1,numEpochs);
epochTime0 = zeros(1,numEpochs);
epochLabels = cell(1,numEpochs);
% Go through each epoch and pull out its begin time and end time, and
% convert from nanoseconds to samples. Calculate number or samples in
% epoch. Get the first and last block for the epoch.
for p = 0:numEpochs-1
    anEpoch = epochListArray.get(p);
    epochBegin = uint64(anEpoch.getBeginTime());
    epochEnd = uint64(anEpoch.getEndTime());
    % Note: end time is first sample NOT in epoch. 
    epochBeginSamps(p+1) = mff_micros2Sample(epochBegin, sampRate);
    epochTime0(p+1) = epochBeginSamps(p+1);
    epochNumSamps(p+1) = mff_micros2Sample(epochEnd, sampRate) - epochBeginSamps(p+1);
    epochFirstBlocks(p+1) = anEpoch.getFirstBlock;
    epochLastBlocks(p+1) = anEpoch.getLastBlock;
    epochLabels{p+1} = 'epoch';
%fprintf('epoch diff %f\n', epochEnd - epochBegin)
end

% assumes the following, which should be true: 1-1 mapping between segments
% and epochs, including quantity and begin times.
% todo?: add checks and error cases 
epochType = 'cnt'; % set epoch type to continuous. This gets overwritten below if the data are segmented. 
totalNumSegs = 0;
% If the data is segmented, go through all the segments, find the
% corresponding epoch, and modify the epochLabels and epochTime0 item for
% the epoch. 
categList = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Categories, 'categories.xml', filePath);
if ~isempty(categList) % if there are some categories

    categListArray = categList.getCategories;
    numCategs = categListArray.size;
    for p = 0:numCategs-1
        aCateg = categListArray.get(p);
        categLabel = aCateg.getName;
        segList = aCateg.getSegments;
        numSegs = segList.size;
        totalNumSegs = totalNumSegs + numSegs;
        for q = 0:numSegs-1
            aSeg = segList.get(q);
            segBegin = uint64(aSeg.getBeginTime);
%segEnd = aSeg.getEndTime;
%fprintf('seg diff %f\n', segEnd - segBegin)
            segBeginSamp = mff_micros2Sample(segBegin, sampRate);
            segInd = find(epochBeginSamps == segBeginSamp);
            epochLabels{segInd} = char(categLabel);
            time0 = uint64(aSeg.getEventBegin());
            time0Samp = mff_micros2Sample(time0, sampRate);
            time0Samp = (time0Samp - segBeginSamp) + 1;
            epochTime0(segInd) = time0Samp;
        end
    end
    epochType = 'seg'; % set the epoch type to segmented. 
    % if epoch lengths are different, yet there are categories, than it's
    % var(iable) epoch type. 
    if size(unique(epochNumSamps),2) ~= 1 || size(unique(epochTime0),2) ~= 1
        epochType = 'var'; % set the epoch type to variable. 
    end
end
%fprintf('totalNumSegs, numEpochs %d %d\n', totalNumSegs, numEpochs);
