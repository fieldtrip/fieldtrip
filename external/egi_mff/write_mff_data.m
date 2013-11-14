%% write_mff_data.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012
%  Copyright 2012 EGI. All rights reserved.
%  Writes eeg data in newData to MFF file.
%
%  This function has several limitations: 
%  -- It doesn't necessarily fit into the Field Trip framework, although it
%  the newData and hdr parameters are Field Trip style variables.  
%  -- It assumes that the data is in the same structure as the srcFile in
%  terms of number of channels, number of samples, number of epochs and
%  sizes of epochs. 
%  -- It doesn't modify the MFF file's history data. 
%
%  If dstFilePath is empty, or the same as srcFilePath, then the data gets
%  overwritten, otherwise, the srcFilePath is copied to dstFilePath, and the
%  data is written to dstFilePath. In the later case, hdr is ignored. 
%%
function write_mff_data(srcFilePath, dstFilePath, newData, hdr)
if isempty(dstFilePath)
    dstFilePath = srcFilePath;
end
if strcmp(dstFilePath, srcFilePath)
    if isempty(hdr)
        srcSummaryInfo = mff_getSummaryInfo(dstFilePath);
    else
        srcSummaryInfo = hdr.orig;
    end
else
    copyfile(srcFilePath, dstFilePath);
    srcSummaryInfo = mff_getSummaryInfo(dstFilePath);
end
if ~strcmp(class(newData), 'single')
    newData = single(newData);
end

EEGFile = mff_getEEGFilename(srcSummaryInfo.mfffileObj);

dstEEGFile = [EEGFile 'Tmp'];
dstURI = [dstFilePath filesep dstEEGFile];
delegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
factory = javaObject('com.egi.services.mff.api.MFFFactory', delegate);
resourceType = javaObject('com.egi.services.mff.api.MFFResourceType', com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal);
factory.createResourceAtURI(dstURI, resourceType);
dstBinObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal, dstEEGFile, dstFilePath);

% dstURI = [dstFilePath filesep EEGFile 'Tmp'];
% dstBinObj = javaObject('com.egi.services.mff.api.Signal');
% dstBinObj.openResource(dstURI);

for p=0:srcSummaryInfo.blocks.size - 1
    aBlock = srcSummaryInfo.blocks.get(p);
%     aBlock = srcSummaryInfo.binObj.loadSignalBlockData(aBlock);
%     aBlock.data = [];
%     java.lang.Runtime.getRuntime.freeMemory;
    
    numChannels = aBlock.numberOfSignals;
    % number of 4 byte floats is 1/4 the data block size
    % That is divided by channel count to get data for each channel:
    samplesTimesChannels = aBlock.dataBlockSize/4;
    numSamples = samplesTimesChannels / numChannels;

    if p == 0
        beginSamp = 1;
    else
        beginSamp = endSamp + 1;
    end
    endSamp = (beginSamp + numSamples)-1;
%     if p > 0
%         beginSamp = sum(srcSummaryInfo.epochNumSamps(1:p)) + 1;
%     end
% %    beginSamp = srcSummaryInfo.epochBeginSamps(p+1)+1;
%     endSamp = (beginSamp + srcSummaryInfo.epochNumSamps(p+1))-1;
    
    
%fprintf('%d %d %d\n', beginSamp, endSamp, size(newData,2));
    newDataEpoch = newData(:,beginSamp:endSamp);
    newDataEpoch = reshape(newDataEpoch', size(newDataEpoch,1) * size(newDataEpoch,2), 1);
    newDataEpoch = typecast(newDataEpoch, 'int8');
    aBlock.data = newDataEpoch;
    dstBinObj.writeSignalBlock(aBlock);
    aBlock.data = [];
    java.lang.Runtime.getRuntime.freeMemory;
end
factory.closeResource(dstBinObj);
%dstBinObj.closeResource();

% delete the EEG file
% rename the written file to the EEG filename
movefile(dstURI, [dstFilePath filesep EEGFile]);
