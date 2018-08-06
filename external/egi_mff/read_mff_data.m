%% read_mff_data.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012, 4/15/2014
%  Copyright 2012, 2014 EGI. All rights reserved.
% 
%  Takes the path to the data and returns the 2-D matrix of size
%  Nchans*Nsamples as described at
%  http://www.fieldtriptoolbox.org/reference/ft_read_data.
% 
%  filePath ? The path to the .mff file. 
%
%  indType ? Indicates how to interpret the next two parameters. There are
%  two possible values:
% 
%  -- 'epoch', which interprets the following two parameters as epoch
%  numbers. 
%  -- 'sample', which interprets the following two parameters as
%  sample numbers.
% 
%  chanInds  ? Indicates which channels to include. If you pass in []
%  (MATLAB?s equivalent of NULL), you get all the channels. You can pass in
%  lists of channels in MATLAB format, such as [1 2 3] or [1:3], and so on.
% 
%  hdr ? FieldTrip header. You have the option of passing in the header, or
%  []. If you pass in the header, it pulls data out of it, rather than
%  recomputing them.
%%
function data = read_mff_data(filePath, indType, beginInd, endInd, chanInds, hdr)
if isempty(hdr)
%     try
        summaryInfo = mff_getSummaryInfo(filePath);
%     catch theException
%         throw(theException)
%     end
else
    summaryInfo = hdr.orig;
end

% If the begin and end are specified in epochs, then need the begin and end
% blocks corresponding to the epochs. Otherwise, it's specified in samples,
% in which case, we also need the begin and end samples. 
if strcmp(indType, 'sample')
    [beginBlock beginSample] = blockSample2BlockAndSample(beginInd, summaryInfo.blockNumSamps);
    [endBlock endSample] = blockSample2BlockAndSample(endInd, summaryInfo.blockNumSamps);
else
    beginBlock = summaryInfo.epochFirstBlocks(beginInd);
    endBlock = summaryInfo.epochLastBlocks(endInd);
end
dataNumSamples = (summaryInfo.blockBeginSamps(endBlock) - summaryInfo.blockBeginSamps(beginBlock)) + summaryInfo.blockNumSamps(endBlock);
pibNChans = 0;
if ~isempty(summaryInfo.javaObjs.pibBinObj)
    pibNChans = summaryInfo.pibNChans;
    if summaryInfo.pibHasRef
        pibNChans = pibNChans + 1;
    end
end
eegNChans = summaryInfo.nChans;
data = zeros(eegNChans + pibNChans, dataNumSamples);
% Get the data from the blocks. 
% EEG data...
% Call to function replaced by inline code for speed purposes. 
% data(1:eegNChans,:) = read_mff_data_blocks(summaryInfo.binObj, summaryInfo.blocks, beginBlock, endBlock, eegNChans, dataNumSamples, summaryInfo.blockNumSamps);
%%
binObj = summaryInfo.javaObjs.binObj;
blocks = summaryInfo.javaObjs.blocks;
startChan = 1;
endChan = eegNChans;
sampleInd = 1;
for blockInd = beginBlock-1:endBlock-1
%     fprintf('blockInd %d\n', blockInd); %!!!
    lastSampleInd = sampleInd + summaryInfo.blockNumSamps(blockInd+1) - 1;
    data(startChan:endChan,sampleInd:lastSampleInd) = read_mff_data_block(binObj, blocks, blockInd);
    sampleInd = lastSampleInd + 1;
end

% PIB data if any...
if ~isempty(summaryInfo.javaObjs.pibBinObj)
% Call to function replaced by inline code for speed purposes. 
%     data(eegNChans+1:end,:) = read_mff_data_blocks(summaryInfo.pibBinObj, summaryInfo.pibBlocks, beginBlock, endBlock, pibNChans, dataNumSamples, summaryInfo.blockNumSamps);
    binObj = summaryInfo.javaObjs.pibBinObj;
    blocks = summaryInfo.javaObjs.pibBlocks;
    startChan = eegNChans + 1;
    endChan = eegNChans + pibNChans;
    sampleInd = 1;
    for blockInd = beginBlock-1:endBlock-1
%         fprintf('blockInd %d\n', blockInd); %!!!
        lastSampleInd = sampleInd + summaryInfo.blockNumSamps(blockInd+1) - 1;
        data(startChan:endChan,sampleInd:lastSampleInd) = read_mff_data_block(binObj, blocks, blockInd);
        sampleInd = lastSampleInd + 1;
    end
end
% if channel indeces were provided, downsample to the requested channels
if size(chanInds,1) ~= 0
    data = data(chanInds,:);
end
% If begin and end are specified in samples, trim the data down to the
% specified samples. 
if strcmp(indType, 'sample')
    if (beginSample ~= 1) || (beginSample + (endInd-beginInd) ~= size(data,2))
        data = data(:,beginSample:beginSample + (endInd-beginInd));
    end
% Otherwise, if the data are segmented, reshape into trials. 
elseif strcmp(summaryInfo.epochType, 'seg')
    nChans = size(data, 1);
    nSamples = summaryInfo.epochNumSamps(1);
    nTrials = (endInd - beginInd) + 1;
    data = reshape(data,nChans, nSamples, nTrials);
end

%% Gets the data from the blocks. 
% Three versions of this routine, from fastest to slowest (newest to
% oldest), preserved for informational purposes. Ultimately not used
% because it's much faster to insert these lines inline in the caller. 
% function data = read_mff_data_blocks(binObj, blocks, beginBlock, endBlock, numChans, numSamples, blockNumSamps)
% data = zeros(numChans, numSamples);
% sampleInd = 1;
% for blockInd = beginBlock-1:endBlock-1
% %     fprintf('blockInd %d\n', blockInd); %!!!
%     lastSampleInd = sampleInd + blockNumSamps(blockInd+1) - 1;
%     data(:,sampleInd:lastSampleInd) = read_mff_data_block(binObj, blocks, blockInd);
%     sampleInd = lastSampleInd + 1;
% end
% function data = read_mff_data_blocks(binObj, blocks, beginBlock, endBlock, numChans, numSamples, blockNumSamps)
% for blockInd = beginBlock-1:endBlock-1
% %     fprintf('blockInd %d\n', blockInd); %!!!
%     tmpdata = read_mff_data_block(binObj, blocks, blockInd);
%     lastSampleInd = sampleInd + size(tmpdata,2) - 1;
%     data(:,sampleInd:lastSampleInd) = tmpdata;
%     sampleInd = lastSampleInd + 1;
% end
% function data = read_mff_data_blocks(binObj, blocks, beginBlock, endBlock, numChans, numSamples, blockNumSamps)
% for blockInd = beginBlock-1:endBlock-1
% %     fprintf('blockInd %d\n', blockInd);
%     tmpdata = read_mff_data_block(binObj, blocks, blockInd);
%     if blockInd == beginBlock-1
%         data = tmpdata;
%     else
%         if size(data,1) == size(tmpdata,1)
%             data = [data tmpdata];
%         else
%             % Error: blocks disagree on number of channels. Should never
%             % occur, especially given the checking performed by this point.
%             % todo?: Add error handling?
%         end
%     end
% end

%% Gets one block of data.  
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

%% Given a sample number within the complete data, return an epoch number,
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

%% Given a sample number within the complete data, return a block number,
% and a sample number within that block. 
% blockNumSamps contains an element for each block - the number of samples
% in that block. 
function [blockNum sample] = blockSample2BlockAndSample(sampleNum, blockNumSamps)
blockNum = 1;
numSamps = blockNumSamps(blockNum);
while sampleNum > numSamps
    blockNum = blockNum + 1;
    numSamps = numSamps + blockNumSamps(blockNum);
end
numSamps = numSamps - blockNumSamps(blockNum);
sample = sampleNum - numSamps;
