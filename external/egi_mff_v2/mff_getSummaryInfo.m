%% mff_getSummaryInfo.m
%  Matlab File
%  author Colin Davey
%  date 3/2/2012, 12/3/2013
%  Copyright 2012, 2013 EGI. All rights reserved.
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
%  Each epoch has a begin and end time (in microseconds since recording
%  start), and first and last block.
% 
%  If the data have been segmented, then there is a categories file that has
%  an array of categories. Each category has an array of segments. There is
%  a 1-1 correspondence between segments and epochs. 
% 
%  Each segment has a begin time and end time in microseconds, which can be
%  used to map the segment to the corresponding epoch. Each segment also has
%  a time zero event with its time in microseconds since recording
%  start. 
%%
function summaryInfo = mff_getSummaryInfo(filePath)
try
    mff_valid(filePath);
catch theException
    throw(theException);
end

% create the MFFFile object
mfffileObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_MFFFile, [], filePath);
%mfffileObj = javaObject('com.egi.services.mff.api.MFFFile', filePath, true);

[pibBinObj pibBlocks pibFilename binObj blocks eegFilename] = getSignalBlocks(mfffileObj, filePath);
numblocks = binObj.getNumberOfBlocks();
if numblocks > 0 % Should always be

    summaryInfo.javaObjs.mfffileObj = mfffileObj;
    summaryInfo.javaObjs.binObj = binObj; % bin obj for the EEG data
    summaryInfo.javaObjs.blocks = blocks; % block objects for the EEG bin obj
    summaryInfo.eegFilename = eegFilename;

    blockObj = blocks.get(0); %zero based
    sampRate = double(blockObj.signalFrequency(1)); % 1 based
    nChans = blockObj.numberOfSignals;
    summaryInfo.sampRate = sampRate;
    summaryInfo.nChans = nChans;

    summaryInfo.javaObjs.pibBinObj = pibBinObj; % bin obj for the PIB data
    summaryInfo.javaObjs.pibBlocks = pibBlocks; % block objects for the PIB bin obj
    summaryInfo.pibNChans = 0;
    summaryInfo.pibFilename = '';
    if ~isempty(pibBinObj)
        summaryInfo.pibFilename = pibFilename;
        pnsSetObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_PNSSet, 'pnsSet.xml', filePath);
        pnsSensors = pnsSetObj.getPNSSensors;
        summaryInfo.pibNChans = pnsSensors.size;
        pibBlockObj = pibBlocks.get(0);
%         summaryInfo.pibNChans = pibBlockObj.numberOfSignals;
        summaryInfo.pibHasRef = false;
%         if (summaryInfo.pibNChans == 7) && (pibBlockObj.numberOfSignals == 8)
        if pibBlockObj.numberOfSignals - summaryInfo.pibNChans == 1
           summaryInfo.pibHasRef = true;
        end
    end
    
    [epochType epochBeginSamps epochNumSamps epochFirstBlocks epochLastBlocks epochLabels epochTime0 multiSubj epochSubjects epochFilenames epochSegStatus] = getEpochInfos(filePath, sampRate);
    summaryInfo.epochType = epochType;
    summaryInfo.epochBeginSamps = epochBeginSamps;
    summaryInfo.epochNumSamps = epochNumSamps;
    summaryInfo.epochFirstBlocks = epochFirstBlocks;
    summaryInfo.epochLastBlocks = epochLastBlocks;
    summaryInfo.epochLabels = epochLabels;
    summaryInfo.epochTime0 = epochTime0;
    summaryInfo.multiSubj = multiSubj;
    summaryInfo.epochSubjects = epochSubjects;
    summaryInfo.epochFilenames = epochFilenames;
    summaryInfo.epochSegStatus = epochSegStatus;

    % This allows data reading by the block, rather than by the entire epoch. 
    summaryInfo.blockBeginSamps = zeros(1,numblocks);
    summaryInfo.blockNumSamps = zeros(1,numblocks);
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
        
        % This chunk of code allows data reading by the block, rather than by the entire epoch. 
        % number of 4 byte floats is 1/4 the data block size
        % That is divided by channel count to get data for each channel:
        sampsTimesChans = blockObj.dataBlockSize/4;
        nSamps = sampsTimesChans / nChans;
%         summaryInfo.blockBeginSamps(x+1) = ;
        summaryInfo.blockNumSamps(x+1) = nSamps;
        if x ~= 0
            summaryInfo.blockBeginSamps(x+1) = summaryInfo.blockBeginSamps(x) + summaryInfo.blockNumSamps(x);
        end
    end
else
    % Error: Signal has 0 blocks. Should never occur. todo?: error
    % handling
end

% Returns the bin objects corresponding to the EEG, and PIB if it exists,
% and the blocks arrays associated associated with them. 
function [pibBinObj pibBlocks pibSignalFile binObj blocks eegFile] = getSignalBlocks(mfffileObj, filePath)
% eegFile = mff_getEEGFilename(mfffileObj, filePath);
eegFile = mff_getSignalFilename(mfffileObj, filePath, com.egi.services.mff.api.InfoN.kEEG);
binObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal, eegFile, filePath);
blocks = binObj.getSignalBlocks();

% pibSignalFile = mff_getPIBFilename(mfffileObj, filePath);
pibSignalFile = mff_getSignalFilename(mfffileObj, filePath, com.egi.services.mff.api.InfoN.kPNSData);
pibBinObj = [];
pibBlocks = [];
if ~isempty(pibSignalFile)
    pibBinObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal, pibSignalFile, filePath);
    pibBlocks = pibBinObj.getSignalBlocks();
end

%  Gets the filename of the bin file with the EEG in it. 
function signalFile = mff_getEEGFilename(mfffileObj, filePath)
signalFile = mff_getSignalFilename(mfffileObj, filePath, com.egi.services.mff.api.InfoN.kEEG);

%  Gets the filename of the bin file with the pib data in it. 
function signalFile = mff_getPIBFilename(mfffileObj, filePath)
signalFile = mff_getSignalFilename(mfffileObj, filePath, com.egi.services.mff.api.InfoN.kPNSData);

function signalFile = mff_getSignalFilename(mfffileObj, filePath, infoNType)
signalFile = [];
%infoFile = [];
binfiles = mfffileObj.getSignalResourceList(false);
% Java is zero based. 
for p = 0:size(binfiles)-1
    binFilename = binfiles.get(p);
    % All this to strip the number (binNumStr) from the signal file in order to apply
    % it to the info file. 
    prefix = 'signal';
    prefixLen = size(prefix,2);
    extension = '.bin';
    extensionLen = size(extension,2);
    binNumDigits = size(binFilename,2) - (prefixLen + extensionLen);
    binNumStr = binFilename(prefixLen+1:prefixLen + binNumDigits);

    infoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_InfoN, ['info' binNumStr '.xml'], filePath);
    fileType = infoObj.getInfoNFileType;

    if fileType == infoNType
%         infoFile = ['info' binNumStr '.xml'];
        signalFile = ['signal' binNumStr '.bin'];
        % can break here - assume only one item. 
        break;
    end
end

% %  Gets the filename of the bin file with the EEG in it. 
% function EEGFile = mff_getEEGFilename(mfffileObj)
% binfiles = mfffileObj.getSignalResourceList(false);
% EEGBinInd = 0; 
% %EEGFile = binfiles.elementAt(EEGBinInd);
% EEGFile = binfiles.get(EEGBinInd);
% 
% %  Gets the filename of the bin file with the pib data in it. 
% function [signalFile infoFile] = mff_getPIBFilenames(mfffileObj, filePath)
% signalFile = [];
% infoFile = [];
% binfiles = mfffileObj.getSignalResourceList(false);
% % Java are zero based, but we are skipping the first/0 one because that's
% % assumed to be EEG. 
% for p = 1:size(binfiles)-1
%     binFilename = binfiles.get(p);
%     % All this to strip the number (binNumStr) from the signal file in order to apply
%     % it to the info file. 
%     prefix = 'signal';
%     prefixLen = size(prefix,2);
%     extension = '.bin';
%     extensionLen = size(extension,2);
%     binNumDigits = size(binFilename,2) - (prefixLen + extensionLen);
%     binNumStr = binFilename(prefixLen+1:prefixLen + binNumDigits);
% 
%     infoObj = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_InfoN, ['info' binNumStr '.xml'], filePath);
%     fileType = infoObj.getInfoNFileType;
% 
%     if fileType == 3 % PNSData
%         infoFile = ['info' binNumStr '.xml'];
%         signalFile = ['signal' binNumStr '.bin'];
%         
% %         fileTypeInformation = infoObj.getInfoNFileTypeInformation;
% %         fileTypeInformation.getPNSSetName
%         % can break here - assume only one pib signal file. 
%         break;
%     end
% end

% Each return value is an array with one element for each epoch, except epochType, which is global description.
function [epochType epochBeginSamps epochNumSamps epochFirstBlocks epochLastBlocks epochLabels epochTime0 multiSubj epochSubjects epochFilenames epochSegStatus] = getEpochInfos(filePath, sampRate)
epochList = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Epochs, 'epochs.xml', filePath);
epochListArray = epochList.getEpochs;
numEpochs = epochListArray.size;
epochBeginSamps = zeros(1,numEpochs);
epochNumSamps = zeros(1,numEpochs);
epochFirstBlocks = zeros(1,numEpochs);
epochLastBlocks = zeros(1,numEpochs);
epochTime0 = zeros(1,numEpochs);
epochLabels = cell(1,numEpochs);
epochSubjects = [];
epochFilenames = [];
epochSegStatus = [];
multiSubj = false;
% Go through each epoch and pull out its begin time and end time, and
% convert from microseconds to samples. Calculate number of samples in
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
% If the data are segmented, go through all the segments, find the
% corresponding epoch, and modify the epochLabels and epochTime0 item for
% the epoch. 
categList = mff_getObject(com.egi.services.mff.api.MFFResourceType.kMFF_RT_Categories, 'categories.xml', filePath);
if ~isempty(categList) % if there are some categories
    epochSegStatus = cell(1,numEpochs);
    categListArray = categList.getCategories;
    numCategs = categListArray.size;
    % The following approximately 30 lines are just to determine if the data
    % are multi subject or not. 
    % To be considered a multi-subject file, it must meet the following
    % conditions: 
    % -- It must contain a subjects directory
    % -- It must have more than one XML file because there is one XML file
    % for each subject. Exception: If it is a result of the Combine tool, then
    % the subjects folder might have zero files. In that case, the subject
    % info is in the keys. 
    % -- No segment may have more than one subject. Otherwise, it is a
    % grand average file, not a multi-subject file. 
%     fprintf('*** before multiSubj test\n');
    if exist([filePath filesep 'subjects']) == 7
        subjectFiles = dir([filePath filesep 'subjects/*.xml']);
        % If there are zero subject files, then assume the subject info is
        % in the keys. 
%         if size(subjectFiles,1) > 1 || size(subjectFiles,1) == 0
            multiSubj = true;
            for p = 0:numCategs-1
                aCateg = categListArray.get(p);
                segList = aCateg.getSegments;
                numSegs = segList.size;
                for q = 0:numSegs-1
                    aSeg = segList.get(q);
                    keyListArray = aSeg.getKeys;
                    numKeys = keyListArray.size;
                    numSubjs = 0;
                    for r = 0:numKeys-1
                        aKey = keyListArray.get(r);
%                         fprintf('%d: %s %s %s %s\n', r, char(aKey.getCode), char(aKey.getData), char(aKey.getDataType), char(aKey.getDescription));
                        if strcmp(char(aKey.getCode), 'subj')
                            numSubjs = numSubjs + 1;
                            if numSubjs == 2
                                multiSubj = false;
                                break;
                            end
                        end
                    end
                    % Because the above break only breaks the inner loop...
                    if ~multiSubj
                        break;
                    end
                end
            end
%         end
    end
    if multiSubj
        epochSubjects = cell(1,numEpochs);
        epochFilenames = cell(1,numEpochs);
    end
%%    
%     fprintf('*** before category loop\n');
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
            epochSegStatus{segInd} = char(aSeg.getStatus);
            epochLabels{segInd} = char(categLabel);
            time0 = uint64(aSeg.getEventBegin());
            time0Samp = mff_micros2Sample(time0, sampRate);
            time0Samp = (time0Samp - segBeginSamp) + 1;
            if time0Samp < 1
                time0Samp = 1;
            end
            epochTime0(segInd) = time0Samp;
            if multiSubj
                keyListArray = aSeg.getKeys;
                numKeys = keyListArray.size;
                for r = 0:numKeys-1
                    aKey = keyListArray.get(r);
%                     fprintf('%d: %s %s %s %s\n', r, char(aKey.getCode), char(aKey.getData), char(aKey.getDataType), char(aKey.getDescription));
                    switch char(aKey.getCode)
                        case 'subj'
                            subject = char(aKey.getData);
                        case 'FILE'
                            filename = char(aKey.getData);
                    end
                end
                epochSubjects{segInd} = subject;
                epochFilenames{segInd} = filename;
            end
        end
    end
    if multiSubj && size(unique(epochSubjects),2) == 1 && size(unique(epochFilenames),2) == 1
        epochSubjects = [];
        epochFilenames = [];
        multiSubj = false;
    end
    epochType = 'seg'; % set the epoch type to segmented. 
    % if epoch lengths are different, yet there are categories, than it's
    % var(iable) epoch type. 
    if size(unique(epochNumSamps),2) ~= 1 || size(unique(epochTime0),2) ~= 1
        epochType = 'var'; % set the epoch type to variable. 
    end
end
%fprintf('totalNumSegs, numEpochs %d %d\n', totalNumSegs, numEpochs);
