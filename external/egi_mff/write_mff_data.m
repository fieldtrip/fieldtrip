%% write_mff_data.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012, 4/15/2014
%  Copyright 2012, 2014 EGI. All rights reserved.
% 
%  Writes channel data (newData) in newData to an MFF file. NewData is a 2-D
%  matrix of size Nchans*Nsamples as described at
%  http://fieldtrip.fcdonders.nl/reference/ft_read_data.
% 
%  If dstMFFPath is empty, or the same as srcMFFPath, then the data gets
%  overwritten, otherwise, the srcMFFPath is copied to dstMFFPath, and the
%  data is written to dstMFFPath. In the later case, hdr is ignored.
% 
%  hdr ? FieldTrip header. You have the option of passing in the header, or
%  []. If you pass in the header, it pulls data out of it, rather than
%  recomputing them.
% 
%  This function has the following limitations: It requires a source data
%  file (srcMFFPath), in other words, it doesn?t allow you to create
%  synthetic data and create a new MFF file. Also, it assumes that the data
%  is in the same structure as the srcFile in terms of number of channels,
%  number of samples, number of epochs and sizes of epochs. So, it doesn?t
%  support processes that change any of those items, for example, channel
%  downsampling, or segmentation.
% 
%  This function doesn't modify the MFF file's history data. To do that,
%  call mff_write_history.
%%
function write_mff_data(srcMFFPath, dstMFFPath, newData, hdr)
if isempty(dstMFFPath)
    dstMFFPath = srcMFFPath;
end
if strcmp(dstMFFPath, srcMFFPath)
    if isempty(hdr)
        srcSummaryInfo = mff_getSummaryInfo(dstMFFPath);
    else
        srcSummaryInfo = hdr.orig;
    end
else
    if exist(dstMFFPath, 'dir') == 7
        fullDS_StorePath = [dstMFFPath filesep '._.DS_Store'];
        if exist(fullDS_StorePath, 'file') == 2
            delete(fullDS_StorePath);
        end
        rmdir(dstMFFPath, 's');
    end
    copyfile(srcMFFPath, dstMFFPath);
    srcSummaryInfo = mff_getSummaryInfo(dstMFFPath);
end
if ~strcmp(class(newData), 'single')
    newData = single(newData);
end
write_mff_signal(srcSummaryInfo.eegFilename, srcSummaryInfo.javaObjs.blocks, dstMFFPath, newData(1:srcSummaryInfo.nChans,:));
if srcSummaryInfo.pibNChans ~= 0
    pibData = newData(srcSummaryInfo.nChans+1:srcSummaryInfo.nChans+srcSummaryInfo.pibNChans,:);
    if srcSummaryInfo.pibHasRef
        numSamples = size(newData, 2);
        pibData = [pibData ; zeros(1, numSamples)];
    end
    write_mff_signal(srcSummaryInfo.pibFilename, srcSummaryInfo.javaObjs.pibBlocks, dstMFFPath, pibData);
end

function write_mff_signal(signalFilename, blocks, dstMFFPath, newData)
dstSignalFile = [signalFilename 'Tmp'];
dstURI = [dstMFFPath filesep dstSignalFile];
delegate = javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate');
factory = javaObject('com.egi.services.mff.api.MFFFactory', delegate);
resourceType = javaObject('com.egi.services.mff.api.MFFResourceType', com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal);
factory.createResourceAtURI(dstURI, resourceType);
dstBinObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal, dstSignalFile, dstMFFPath);

for p=0:blocks.size - 1
    aBlock = blocks.get(p);
    
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
    
    
% fprintf('%d: %d %d %d\n', p, beginSamp, endSamp, size(newData,2));
    newDataEpoch = newData(:,beginSamp:endSamp);
    newDataEpoch = reshape(newDataEpoch', size(newDataEpoch,1) * size(newDataEpoch,2), 1);
    newDataEpoch = typecast(newDataEpoch, 'int8');
    aBlock.data = newDataEpoch;
    dstBinObj.writeSignalBlock(aBlock);
    aBlock.data = [];
    java.lang.Runtime.getRuntime.freeMemory;
end
factory.closeResource(dstBinObj);

% delete the EEG file
% rename the written file to the EEG filename
movefile(dstURI, [dstMFFPath filesep signalFilename]);
