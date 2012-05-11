%% read_mff_data.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012
%  Copyright 2012 EGI. All rights reserved.
%  Returns the data in Field Trip format. 
%%
function data = read_mff_data(filePath, indType, beginInd, endInd, chanInds, hdr)
if isempty(hdr)
    summaryInfo = mff_getSummaryInfo(filePath);
else
    summaryInfo = hdr.orig;
end

% If the begin and end are specified in epochs, then need the begin and end
% blocks corresponding to the epochs. Otherwise, it's specified in samples,
% in which case, we also need the begin and end samples. 
if strcmp(indType, 'sample')
    [beginEpoch beginSample] = epochSample2EpochAndSample(beginInd, summaryInfo.epochNumSamps);
    [endEpoch endSample] = epochSample2EpochAndSample(endInd, summaryInfo.epochNumSamps);
    beginBlock = summaryInfo.epochFirstBlocks(beginEpoch);
    endBlock = summaryInfo.epochLastBlocks(endEpoch);
else
    beginBlock = summaryInfo.epochFirstBlocks(beginInd);
    endBlock = summaryInfo.epochLastBlocks(endInd);
end
% Get the data from the blocks. 
data = read_mff_data_blocks(summaryInfo.binObj, summaryInfo.blocks, beginBlock, endBlock);
% if channel indeces were provided, downsample to the requested channels
if size(chanInds,1) ~= 0
    data = data(chanInds,:);
end
% If begin and end are specified in samples, trim the data down to the
% specified samples. 
if strcmp(indType, 'sample')
    data = data(:,beginSample:beginSample + (endInd-beginInd));
% Otherwise, if the data are segmented, reshape into trials. 
elseif strcmp(summaryInfo.epochType, 'seg')
    nChans = size(data, 1);
    nSamples = summaryInfo.epochNumSamps(1);
    nTrials = (endInd - beginInd) + 1;
    data = reshape(data,nChans, nSamples, nTrials);
end

% Gets the data from the blocks. 
function data = read_mff_data_blocks(binObj, blocks, beginBlock, endBlock)
for blockInd = beginBlock-1:endBlock-1
    tmpdata = read_mff_data_block(binObj, blocks, blockInd);
    if blockInd == beginBlock-1
        data = tmpdata;
    else
        if size(data,1) == size(tmpdata,1)
            data = [data tmpdata];
        else
            % Error: blocks disagree on number of channels. Should never
            % occur, especially given the checking performed by this point.
            % todo?: Add error handling?
        end
    end
end

% Gets one block of data.  
function data = read_mff_data_block(binObj, blocks, blockInd)
% Load a block
blockObj = blocks.get(blockInd);
blockObj = binObj.loadSignalBlockData(blockObj);

% The data is stored as bytes. Need numbers for reshaping. 
numChannels = blockObj.numberOfSignals;
% number of 4 byte floats is 1/4 the data block size
% That is divided by channel count to get data for each channel:
samplesTimesChannels = blockObj.dataBlockSize/4;
numSamples = samplesTimesChannels / numChannels;

% get block, returned as bytes. 
data = blockObj.data;
% free up memory in java space. 
blockObj.data = [];
java.lang.Runtime.getRuntime.freeMemory;
% convert bytes to equivalent floating point values
data = typecast(data,'single');
% reshape as an epoch. 
data = reshape(data, numSamples, numChannels)';

% Given a sample number within the complete data, return an epoch number,
% and a sample number within that epoch. 
% epochNumSamps contains an element for each epoch - the number of samples
% in that epoch. 
function [epochNum sample] = epochSample2EpochAndSample(sampleNum, epochNumSamps)
epochNum = 1;
numSamps = epochNumSamps(epochNum);
while sampleNum > numSamps
    epochNum = epochNum + 1;
    numSamps = numSamps + epochNumSamps(epochNum);
end
numSamps = numSamps - epochNumSamps(epochNum);
sample = sampleNum - numSamps;
