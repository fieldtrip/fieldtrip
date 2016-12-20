function output = read_nervus_header(filename)
% read_nervus_header  Returns header information from Nicolet file.
%
%   FILENAME is the file name of a file in the Natus/Nicolet/Nervus(TM)
%   format (originally designed by Taugagreining HF in Iceland)
%
%   Based on ieeg-portal/Nicolet-Reader
%   at https://github.com/ieeg-portal/Nicolet-Reader
%
% Copyright (C) 2016, Jan Brogger and Joost Wagenaar 
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: $

%--Constants--
LABELSIZE = 32;
TSLABELSIZE = 64;
UNITSIZE = 16;
ITEMNAMESIZE  = 64;

% ---------------- Opening File------------------
h = fopen(filename,'rb','ieee-le');
if h==-1
    error('Can''t open Nervus EEG file')
end

nrvHdr = struct();
nrvHdr.filename = filename;
nrvHdr.misc1 = fread(h,5, 'uint32');
nrvHdr.unknown = fread(h,1,'uint32');
nrvHdr.indexIdx = fread(h,1,'uint32');
[nrvHdr.NrStaticPackets, nrvHdr.StaticPackets] = read_nervus_header_staticpackets(h);
nrvHdr.QIIndex = read_nervus_header_Qi(h, nrvHdr.NrStaticPackets);
nrvHdr.QIIndex2 = read_nervus_header_Qi2(h, nrvHdr.QIIndex);
nrvHdr.MainIndex  = read_nervus_header_main(h, nrvHdr.indexIdx, nrvHdr.QIIndex.nrEntries);
nrvHdr.allIndexIDs = [nrvHdr.MainIndex.sectionIdx];
nrvHdr.infoGuids = read_nervus_header_infoGuids(h, nrvHdr.StaticPackets, nrvHdr.MainIndex);
nrvHdr.DynamicPackets = read_nervus_header_dynamicpackets(h, nrvHdr.StaticPackets, nrvHdr.MainIndex);
nrvHdr.PatientInfo = read_nervus_header_patient(h, nrvHdr.StaticPackets, nrvHdr.MainIndex);
nrvHdr.SigInfo = read_nervus_header_SignalInfo(h, nrvHdr.StaticPackets, nrvHdr.MainIndex, ITEMNAMESIZE, LABELSIZE, UNITSIZE);
nrvHdr.ChannelInfo = read_nervus_header_ChannelInfo(h, nrvHdr.StaticPackets, nrvHdr.MainIndex, ITEMNAMESIZE, LABELSIZE);
nrvHdr.TSInfo = read_nervus_header_TSInfo(nrvHdr.DynamicPackets, TSLABELSIZE, LABELSIZE);
nrvHdr.Segments = read_nervus_header_Segments(h, nrvHdr.StaticPackets, nrvHdr.MainIndex, nrvHdr.TSInfo);
nrvHdr.Events = read_nervus_header_events(h, nrvHdr.StaticPackets, nrvHdr.MainIndex);
nrvHdr.MontageInfo = read_nervus_header_montage(h, nrvHdr.StaticPackets, nrvHdr.MainIndex);
nrvHdr.MontageInfo2 = read_nervus_header_dynamic_montages(nrvHdr.DynamicPackets);

reference = unique(nrvHdr.Segments(1).refName(cellfun(@length, [nrvHdr.Segments(1).refName])>0));
if strcmp(reference, 'REF')
	nrvHdr.reference = 'common';
else
	nrvHdr.reference = 'unknown';
end
fclose(h);

%Calculate sample count across segments
% - some channels have lower sampling rates, so we for each segments we
%   choose the channel with the highest sampling rate
totalNSamples = 0;
for i=1:size(nrvHdr.Segments,2)
    totalNSamples = totalNSamples + max(nrvHdr.Segments(i).samplingRate*nrvHdr.Segments(i).duration);
end

output = struct();
output.Fs          = max([nrvHdr.Segments.samplingRate]);
output.nChans      = size([nrvHdr.Segments(1).chName],2);
output.label       = nrvHdr.Segments(1).chName;
output.nSamples    = totalNSamples;
output.nSamplesPre = 0;
output.nTrials     = 1; %size(nrvHdr.Segments,2);
output.reference   = nrvHdr.reference;
output.filename    = nrvHdr.filename;
output.orig        = nrvHdr;


end

function [NrStaticPackets, StaticPackets] = read_nervus_header_staticpackets(h)
% Get StaticPackets structure and Channel IDS
fseek(h, 172,'bof');
NrStaticPackets = fread(h,1, 'uint32');
StaticPackets = struct();
for i = 1:NrStaticPackets
    StaticPackets(i).tag = deblank(cast(fread(h, 40, 'uint16'),'char')');
    StaticPackets(i).index = fread(h,1,'uint32');
    switch StaticPackets(i).tag
        case 'ExtraDataStaticPackets'
            StaticPackets(i).IDStr = 'ExtraDataStaticPackets';
        case 'SegmentStream'
            StaticPackets(i).IDStr = 'SegmentStream';
        case 'DataStream'
            StaticPackets(i).IDStr = 'DataStream';
        case 'InfoChangeStream'
            StaticPackets(i).IDStr = 'InfoChangeStream';
        case 'InfoGuids'
            StaticPackets(i).IDStr = 'InfoGuids';
        case '{A271CCCB-515D-4590-B6A1-DC170C8D6EE2}'
            StaticPackets(i).IDStr = 'TSGUID';
        case '{8A19AA48-BEA0-40D5-B89F-667FC578D635}'
            StaticPackets(i).IDStr = 'DERIVATIONGUID';
        case '{F824D60C-995E-4D94-9578-893C755ECB99}'
            StaticPackets(i).IDStr = 'FILTERGUID';
        case '{02950361-35BB-4A22-9F0B-C78AAA5DB094}'
            StaticPackets(i).IDStr = 'DISPLAYGUID';
        case '{8E9421-70F5-11D3-8F72-00105A9AFD56}'
            StaticPackets(i).IDStr = 'FILEINFOGUID';
        case '{E4138BC0-7733-11D3-8685-0050044DAAB1}'
            StaticPackets(i).IDStr = 'SRINFOGUID';
        case '{C728E565-E5A0-4419-93D2-F6CFC69F3B8F}'
            StaticPackets(i).IDStr = 'EVENTTYPEINFOGUID';
        case '{D01B34A0-9DBD-11D3-93D3-00500400C148}'
            StaticPackets(i).IDStr = 'AUDIOINFOGUID';
        case '{BF7C95EF-6C3B-4E70-9E11-779BFFF58EA7}'
            StaticPackets(i).IDStr = 'CHANNELGUID';
        case '{2DEB82A1-D15F-4770-A4A4-CF03815F52DE}'
            StaticPackets(i).IDStr = 'INPUTGUID';
        case '{5B036022-2EDC-465F-86EC-C0A4AB1A7A91}'
            StaticPackets(i).IDStr = 'INPUTSETTINGSGUID';
        case '{99A636F2-51F7-4B9D-9569-C7D45058431A}'
            StaticPackets(i).IDStr = 'PHOTICGUID';
        case '{55C5E044-5541-4594-9E35-5B3004EF7647}'
            StaticPackets(i).IDStr = 'ERRORGUID';
        case '{223A3CA0-B5AC-43FB-B0A8-74CF8752BDBE}'
            StaticPackets(i).IDStr = 'VIDEOGUID';
        case '{0623B545-38BE-4939-B9D0-55F5E241278D}'
            StaticPackets(i).IDStr = 'DETECTIONPARAMSGUID';
        case '{CE06297D-D9D6-4E4B-8EAC-305EA1243EAB}'
            StaticPackets(i).IDStr = 'PAGEGUID';
        case '{782B34E8-8E51-4BB9-9701-3227BB882A23}'
            StaticPackets(i).IDStr = 'ACCINFOGUID';
        case '{3A6E8546-D144-4B55-A2C7-40DF579ED11E}'
            StaticPackets(i).IDStr = 'RECCTRLGUID';
        case '{D046F2B0-5130-41B1-ABD7-38C12B32FAC3}'
            StaticPackets(i).IDStr = 'GUID TRENDINFOGUID';
        case '{CBEBA8E6-1CDA-4509-B6C2-6AC2EA7DB8F8}'
            StaticPackets(i).IDStr = 'HWINFOGUID';
        case '{E11C4CBA-0753-4655-A1E9-2B2309D1545B}'
            StaticPackets(i).IDStr = 'VIDEOSYNCGUID';
        case '{B9344241-7AC1-42B5-BE9B-B7AFA16CBFA5}'
            StaticPackets(i).IDStr = 'SLEEPSCOREINFOGUID';
        case '{15B41C32-0294-440E-ADFF-DD8B61C8B5AE}'
            StaticPackets(i).IDStr = 'FOURIERSETTINGSGUID';
        case '{024FA81F-6A83-43C8-8C82-241A5501F0A1}'
            StaticPackets(i).IDStr = 'SPECTRUMGUID';
        case '{8032E68A-EA3E-42E8-893E-6E93C59ED515}'
            StaticPackets(i).IDStr = 'SIGNALINFOGUID';
        case '{30950D98-C39C-4352-AF3E-CB17D5B93DED}'
            StaticPackets(i).IDStr = 'SENSORINFOGUID';
        case '{F5D39CD3-A340-4172-A1A3-78B2CDBCCB9F}'
            StaticPackets(i).IDStr = 'DERIVEDSIGNALINFOGUID';
        case '{969FBB89-EE8E-4501-AD40-FB5A448BC4F9}'
            StaticPackets(i).IDStr = 'ARTIFACTINFOGUID';
        case '{02948284-17EC-4538-A7FA-8E18BD65E167}'
            StaticPackets(i).IDStr = 'STUDYINFOGUID';
        case '{D0B3FD0B-49D9-4BF0-8929-296DE5A55910}'
            StaticPackets(i).IDStr = 'PATIENTINFOGUID';
        case '{7842FEF5-A686-459D-8196-769FC0AD99B3}'
            StaticPackets(i).IDStr = 'DOCUMENTINFOGUID';
        case '{BCDAEE87-2496-4DF4-B07C-8B4E31E3C495}'
            StaticPackets(i).IDStr = 'USERSINFOGUID';
        case '{B799F680-72A4-11D3-93D3-00500400C148}'
            StaticPackets(i).IDStr = 'EVENTGUID';
        case '{AF2B3281-7FCE-11D2-B2DE-00104B6FC652}'
            StaticPackets(i).IDStr = 'SHORTSAMPLESGUID';
        case '{89A091B3-972E-4DA2-9266-261B186302A9}'
            StaticPackets(i).IDStr = 'DELAYLINESAMPLESGUID';
        case '{291E2381-B3B4-44D1-BB77-8CF5C24420D7}'
            StaticPackets(i).IDStr = 'GENERALSAMPLESGUID';
        case '{5F11C628-FCCC-4FDD-B429-5EC94CB3AFEB}'
            StaticPackets(i).IDStr = 'FILTERSAMPLESGUID';
        case '{728087F8-73E1-44D1-8882-C770976478A2}'
            StaticPackets(i).IDStr = 'DATEXDATAGUID';
        case '{35F356D9-0F1C-4DFE-8286-D3DB3346FD75}'
            StaticPackets(i).IDStr = 'TESTINFOGUID';
            
        otherwise
            if isstrprop(StaticPackets(i).tag, 'digit')
                StaticPackets(i).IDStr = num2str(StaticPackets(i).tag);
            else
                StaticPackets(i).IDStr = 'UNKNOWN';
            end
    end
end
end

function QIIndex = read_nervus_header_Qi(h, nrStaticPackets)
%% QI index
fseek(h, 172208,'bof');
QIIndex =struct();
QIIndex.nrEntries = fread(h,1,'uint32');
QIIndex.misc1 = fread(h,1,'uint32');
QIIndex.indexIdx = fread(h,1,'uint32');
QIIndex.misc3 = fread(h,1,'uint32');
QIIndex.LQi = fread(h,1,'uint64')';
QIIndex.firstIdx = fread(h,nrStaticPackets,'uint64');
end

function QIIndex2 = read_nervus_header_Qi2(h, QIIndex)
fseek(h, 188664,'bof');
QIIndex2  = struct();
for i = 1:QIIndex.LQi
    QIIndex2(i).ftel = ftell(h);
    QIIndex2(i).index = fread(h,2,'uint16')';  %4
    QIIndex2(i).misc1 = fread(h,1,'uint32');   %8
    QIIndex2(i).indexIdx = fread(h,1,'uint32'); %12
    QIIndex2(i).misc2 = fread(h,3,'uint32')'; %24
    QIIndex2(i).sectionIdx = fread(h,1,'uint32');%28
    QIIndex2(i).misc3 = fread(h,1,'uint32'); %32
    QIIndex2(i).offset = fread(h,1,'uint64'); % 40
    QIIndex2(i).blockL = fread(h,1,'uint32');%44
    QIIndex2(i).dataL = fread(h,1,'uint32')';%48
end
end

function MainIndex = read_nervus_header_main(h, indexIdx, nrEntries)
%% Get Main Index:
% Index consists of multiple blocks, after each block is the pointer
% to the next block. Total number of entries is in obj.Qi.nrEntries

MainIndex = struct();
curIdx = 0;
nextIndexPointer = indexIdx;
curIdx2 = 1;
while curIdx < nrEntries
    
    fseek(h, nextIndexPointer, 'bof');
    nrIdx = fread(h,1, 'uint64');
    MainIndex(curIdx + nrIdx).sectionIdx = 0;   % Preallocate next set of indices
    var = fread(h,3*nrIdx, 'uint64');
    for i = 1: nrIdx
        MainIndex(curIdx + i).sectionIdx = var(3*(i-1)+1);
        MainIndex(curIdx + i).offset = var(3*(i-1)+2);
        MainIndex(curIdx + i).blockL = mod(var(3*(i-1)+3),2^32);
        MainIndex(curIdx + i).sectionL = round(var(3*(i-1)+3)/2^32);
    end
    nextIndexPointer = fread(h,1, 'uint64');
    curIdx = curIdx + i;
    curIdx2=curIdx2+1;
end
end

function infoGuids = read_nervus_header_infoGuids(h, StaticPackets, MainIndex)


infoIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'InfoGuids'),1)).index;
indexInstance = MainIndex(find([MainIndex.sectionIdx]==infoIdx,1));
nrInfoGuids = indexInstance.sectionL/16;
infoGuids = struct();
fseek(h, indexInstance.offset,'bof');
for i = 1:nrInfoGuids
    guidmixed = fread(h,16, 'uint8')';
    guidnonmixed = [guidmixed(04), guidmixed(03), guidmixed(02), guidmixed(01), ...
        guidmixed(06), guidmixed(05), guidmixed(08), guidmixed(07), ...
        guidmixed(09), guidmixed(10), guidmixed(11), guidmixed(12), ...
        guidmixed(13), guidmixed(15), guidmixed(15), guidmixed(16)];
    infoGuids(i).guid = num2str(guidnonmixed,'%02X');
end
end


function dynamicPackets = read_nervus_header_dynamicpackets(h, StaticPackets, MainIndex)
dynamicPackets = struct();
indexIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'InfoChangeStream'),1)).index;
offset = MainIndex(indexIdx).offset;
nrDynamicPackets = MainIndex(indexIdx).sectionL / 48;
fseek(h, offset, 'bof');

%Read first only the dynamic packets structure without actual data
for i = 1: nrDynamicPackets
    dynamicPackets(i).offset = offset+i*48;
    guidmixed = fread(h,16, 'uint8')';
    guidnonmixed = [guidmixed(04), guidmixed(03), guidmixed(02), guidmixed(01), ...
        guidmixed(06), guidmixed(05), guidmixed(08), guidmixed(07), ...
        guidmixed(09), guidmixed(10), guidmixed(11), guidmixed(12), ...
        guidmixed(13), guidmixed(14), guidmixed(15), guidmixed(16)];
    dynamicPackets(i).guid = num2str(guidnonmixed, '%02X');
    dynamicPackets(i).guidAsStr = sprintf('{%02X%02X%02X%02X-%02X%02X-%02X%02X-%02X%02X-%02X%02X%02X%02X%02X%02X}', guidnonmixed);
    dynamicPackets(i).date = datenum(1899,12,31) + fread(h,1,'double');
    dynamicPackets(i).datefrac = fread(h,1,'double');
    dynamicPackets(i).internalOffsetStart = fread(h,1, 'uint64')';
    dynamicPackets(i).packetSize = fread(h,1, 'uint64')';
    dynamicPackets(i).data = zeros(0, 1,'uint8');
    
    switch dynamicPackets(i).guid
        case 'BF7C95EF6C3B4E709E11779BFFF58EA7'
            dynamicPackets(i).IDStr = 'CHANNELGUID';
        case '8A19AA48BEA040D5B89F667FC578D635'
            dynamicPackets(i).IDStr = 'DERIVATIONGUID';
        case 'F824D60C995E4D949578893C755ECB99'
            dynamicPackets(i).IDStr = 'FILTERGUID';
        case '0295036135BB4A229F0BC78AAA5DB094'
            dynamicPackets(i).IDStr = 'DISPLAYGUID';
        case '782B34E88E514BB997013227BB882A23'
            dynamicPackets(i).IDStr = 'ACCINFOGUID';
        case 'A271CCCB515D4590B6A1DC170C8D6EE2'
            dynamicPackets(i).IDStr = 'TSGUID';
        case 'D01B34A09DBD11D393D300500400C148'
            dynamicPackets(i).IDStr = 'AUDIOINFOGUID';
        otherwise
            dynamicPackets(i).IDStr = 'UNKNOWN';
    end
end

%Then read the actual data from the pointers above
for i = 1: nrDynamicPackets
    %Look up the GUID of this dynamic packet in the static packets
    % to find the section index
    
    infoIdx = StaticPackets(find(strcmp({StaticPackets.tag},dynamicPackets(i).guidAsStr),1)).index;
    
    %Matching index segments
    indexInstances = MainIndex([MainIndex.sectionIdx] == infoIdx);
    
    %Then, treat all these sections as one contiguous memory block
    % and grab this packet across these instances
    
    internalOffset = 0;
    remainingDataToRead = dynamicPackets(i).packetSize;
    %disp(['Target packet ' dynamicPackets(i).IDStr ' : ' num2str(dynamicPackets(i).internalOffsetStart) ' to ' num2str(dynamicPackets(i).internalOffsetStart+dynamicPackets(i).packetSize) ' target read length ' num2str(remainingDataToRead)]);
    currentTargetStart = dynamicPackets(i).internalOffsetStart;
    for j = 1: size(indexInstances,2)
        currentInstance = indexInstances(j);
        
        %hitInThisSegment = '';
        if (internalOffset <= currentTargetStart) && (internalOffset+currentInstance.sectionL) >= currentTargetStart
            
            startAt = currentTargetStart;
            stopAt =  min(startAt+remainingDataToRead, internalOffset+currentInstance.sectionL);
            readLength = stopAt-startAt;
            
            filePosStart = currentInstance.offset+startAt-internalOffset;
            fseek(h,filePosStart, 'bof');
            dataPart = fread(h,readLength,'uint8=>uint8');
            dynamicPackets(i).data = cat(1, dynamicPackets(i).data, dataPart);
            
            %hitInThisSegment = ['HIT at  ' num2str(startAt) ' to ' num2str(stopAt)];
            %if (readLength < remainingDataToRead)
            %    hitInThisSegment = [hitInThisSegment ' (partial ' num2str(readLength) ' )'];
            %else
            %    hitInThisSegment = [hitInThisSegment ' (finished - this segment contributed ' num2str(readLength) ' )'];
            %end
            %hitInThisSegment = [hitInThisSegment ' abs file pos ' num2str(filePosStart) ' - ' num2str(filePosStart+readLength)];
            
            remainingDataToRead = remainingDataToRead-readLength;
            currentTargetStart = currentTargetStart + readLength;
            
        end
        %disp(['    Index ' num2str(j) ' Offset: ' num2str(internalOffset) ' to ' num2str(internalOffset+currentInstance.sectionL) ' ' num2str(hitInThisSegment)]);
        
        internalOffset = internalOffset + currentInstance.sectionL;
    end
end
end

function PatientInfo = read_nervus_header_patient(h, StaticPackets, Index)
%% Get PatientGUID
PatientInfo = struct();

infoProps = { 'patientID', 'firstName','middleName','lastName',...
    'altID','mothersMaidenName','DOB','DOD','street','sexID','phone',...
    'notes','dominance','siteID','suffix','prefix','degree','apartment',...
    'city','state','country','language','height','weight','race','religion',...
    'maritalStatus'};

infoIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'PATIENTINFOGUID'),1)).index;
indexInstance = Index(find([Index.sectionIdx]==infoIdx,1));
fseek(h, indexInstance.offset,'bof');
guid = fread(h, 16, 'uint8');
lSection = fread(h, 1, 'uint64');
% reserved = fread(h, 3, 'uint16');
nrValues = fread(h,1,'uint64');
nrBstr = fread(h,1,'uint64');

for i = 1:nrValues
    id = fread(h,1,'uint64');
    switch id
        case {7,8}
            unix_time = (fread(h,1, 'double')*(3600*24)) - 2209161600;% 2208988800; %8
            obj.segments(i).dateStr = datestr(unix_time/86400 + datenum(1970,1,1));
            value = datevec( obj.segments(i).dateStr );
            value = value([3 2 1]);
        case {23,24}
            value = fread(h,1,'double');
        otherwise
            value = 0;
    end
    PatientInfo.(infoProps{id}) = value;
end

strSetup = fread(h,nrBstr*2,'uint64');

for i=1:2:(nrBstr*2)
    id  = strSetup(i);
    value = deblank(cast(fread(h, strSetup(i+1) + 1, 'uint16'),'char')');
    info.(infoProps{id}) = value;
end

end

function sigInfo = read_nervus_header_SignalInfo(h, StaticPackets, Index, ITEMNAMESIZE, LABELSIZE, UNITSIZE)
infoIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'InfoGuids'),1)).index;
indexInstance = Index(find([Index.sectionIdx]==infoIdx,1));
fseek(h, indexInstance.offset,'bof');

sigInfo = struct();
SIG_struct = struct();
sensorIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'SIGNALINFOGUID'),1)).index;
indexInstance = Index(find([Index.sectionIdx]==sensorIdx,1));
fseek(h, indexInstance.offset,'bof');
SIG_struct.guid = fread(h, 16, 'uint8');
SIG_struct.name = fread(h, ITEMNAMESIZE, '*char');
unkown = fread(h, 152, '*char');         %#ok<NASGU>
fseek(h, 512, 'cof');
nrIdx = fread(h,1, 'uint16');  %783
misc1 = fread(h,3, 'uint16'); %#ok<NASGU>

for i = 1: nrIdx
    sigInfo(i).sensorName = deblank(cast(fread(h, LABELSIZE, 'uint16'),'char')');
    sigInfo(i).transducer = deblank(cast(fread(h, UNITSIZE, 'uint16'),'char')');
    sigInfo(i).guid = fread(h, 16, '*uint8');
    sigInfo(i).bBiPolar = logical(fread(h, 1 ,'uint32'));
    sigInfo(i).bAC = logical(fread(h, 1 ,'uint32'));
    sigInfo(i).bHighFilter = logical(fread(h, 1 ,'uint32'));
    sigInfo(i).color =  fread(h, 1 ,'uint32');
    reserved = fread(h, 256, '*char'); %#ok<NASGU>
end
end


function channelInfo = read_nervus_header_ChannelInfo(h, StaticPackets, Index, ITEMNAMESIZE, LABELSIZE)
%% Get CHANNELINFO (CHANNELGUID)
CH_struct = struct();
sensorIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'CHANNELGUID'),1)).index;
indexInstance = Index(find([Index.sectionIdx]==sensorIdx,1));
fseek(h, indexInstance.offset,'bof');
CH_struct.guid = fread(h, 16, 'uint8');
CH_struct.name = fread(h, ITEMNAMESIZE, '*char');
fseek(h, 152, 'cof');
CH_struct.reserved = fread(h, 16, 'uint8');
CH_struct.deviceID = fread(h, 16, 'uint8');
fseek(h, 488, 'cof');

nrIdx = fread(h,2, 'uint32');  %783
channelInfo = struct();
for i = 1: nrIdx(2)
    channelInfo(i).sensor = deblank(cast(fread(h, LABELSIZE, 'uint16'),'char')');
    channelInfo(i).samplingRate = fread(h,1,'double');
    channelInfo(i).bOn = logical(fread(h, 1 ,'uint32'));
    channelInfo(i).lInputID = fread(h, 1 ,'uint32');
    channelInfo(i).lInputSettingID = fread(h,1,'uint32');
    channelInfo(i).reserved = fread(h,4,'char');
    fseek(h, 128, 'cof');
end

curIdx = 0;
for i = 1: length(channelInfo)
    if channelInfo(i).bOn
        channelInfo(i).indexID = curIdx;
        curIdx = curIdx+1;
    else
        channelInfo(i).indexID = -1;
    end
end
end

function [TSInfo] = read_nervus_header_TSInfo(DynamicPackets, TSLABELSIZE, LABELSIZE)
tsPackets = DynamicPackets(strcmp({DynamicPackets.IDStr},'TSGUID'));

if isempty(tsPackets)
    error(['No TSINFO found']);
end    

tsPacket = tsPackets(1);
TSInfo = read_nervus_header_one_TSInfo(tsPacket, TSLABELSIZE, LABELSIZE);

if length(tsPackets) > 1
    allEqual = 1;
    for i = 2: size(tsPackets,2)
        nextTsPacket = tsPackets(i);
        nextTSInfo = read_nervus_header_one_TSInfo(nextTsPacket, TSLABELSIZE, LABELSIZE);       
        areEqual = compareTsInfoPackets(TSInfo, nextTSInfo);
        if (areEqual == 0)
            allEqual = 0;
            break;
        end
    end    
    if (allEqual == 0)            
        error('Multiple TSInfo packets found and they are not the same.');
    end
end
end

function [TSInfo] = read_nervus_header_one_TSInfo(tsPacket, TSLABELSIZE, LABELSIZE)    
    TSInfo = struct();
    elems = typecast(tsPacket.data(753:756),'uint32');
    %alloc = typecast(tsPacket.data(757:760),'uint32');
    
    offset = 761;
    for i = 1:elems
        internalOffset = 0;
        TSInfo(i).label = deblank(char(typecast(tsPacket.data(offset:(offset+TSLABELSIZE-1))','uint16')));
        internalOffset = internalOffset + TSLABELSIZE*2;
        TSInfo(i).activeSensor = deblank(char(typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+LABELSIZE))','uint16')));
        internalOffset = internalOffset + TSLABELSIZE;
        TSInfo(i).refSensor = deblank(char(typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','uint16')));
        internalOffset = internalOffset + 8;
        internalOffset = internalOffset + 56;
        TSInfo(i).lowcut = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        internalOffset = internalOffset + 8;
        TSInfo(i).hiCut = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        internalOffset = internalOffset + 8;
        TSInfo(i).samplingRate = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        internalOffset = internalOffset + 8;
        TSInfo(i).resolution = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        internalOffset = internalOffset + 8;
        TSInfo(i).specialMark = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+2))','uint16');
        internalOffset = internalOffset + 2;
        TSInfo(i).notch = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+2))','uint16');
        internalOffset = internalOffset + 2;
        TSInfo(i).eeg_offset = typecast(tsPacket.data(offset+internalOffset:(offset+internalOffset-1+8))','double');
        offset = offset + 552;
        %disp([num2str(i) ' : ' TSInfo(i).label ' : ' TSInfo(i).activeSensor ' : ' TSInfo(i).refSensor ' : ' num2str(TSInfo(i).samplingRate)]);
    end

end

function areEqual = compareTsInfoPackets(TSInfo1, TSInfo2)    
    areEqual = 1;
    if (size(TSInfo1,2) ~= size(TSInfo2,2))
        areEqual = 0;
    else
        for i = 1:size(TSInfo1,2)
            if (~strcmp(TSInfo1(i).label,TSInfo2(i).label))
                areEqual = 0;
                break;
            end
            if (~strcmp(TSInfo1(i).activeSensor,TSInfo2(i).activeSensor))
                areEqual = 0;
                break;
            end
            if (~strcmp(TSInfo1(i).refSensor,TSInfo2(i).refSensor))
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).lowcut ~= TSInfo2(i).lowcut)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).hiCut ~= TSInfo2(i).hiCut)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).samplingRate ~= TSInfo2(i).samplingRate)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).resolution ~= TSInfo2(i).resolution)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).specialMark ~= TSInfo2(i).specialMark)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).notch ~= TSInfo2(i).notch)
                areEqual = 0;
                break;
            end
            if (TSInfo1(i).eeg_offset ~= TSInfo2(i).eeg_offset)
                areEqual = 0;
                break;
            end
        end
    end
        
end

function [segments] = read_nervus_header_Segments(h, StaticPackets, Index, TSInfo)
%% Get Segment Start Times
segmentIdx = StaticPackets(find(strcmp({StaticPackets.IDStr}, 'SegmentStream'),1)).index;
indexIdx = find([Index.sectionIdx] == segmentIdx, 1);
segmentInstance = Index(indexIdx);

nrSegments = segmentInstance.sectionL/152;
fseek(h, segmentInstance.offset,'bof');
segments = struct();
for i = 1: nrSegments
    dateOLE = fread(h,1, 'double');
    segments(i).dateOLE = dateOLE;
    unix_time = (dateOLE*(3600*24)) - 2209161600;% 2208988800; %8
    segments(i).dateStr = datestr(unix_time/86400 + datenum(1970,1,1));
    datev = datevec( segments(i).dateStr );
    segments(i).startDate = datev(1:3);
    segments(i).startTime = datev(4:6);
    fseek(h, 8 , 'cof'); %16
    segments(i).duration = fread(h,1, 'double');%24
    fseek(h, 128 , 'cof'); %152       
end

% Get nrValues per segment and channel
for iSeg = 1:length(segments)
    % Add Channel Names to segments
    segments(iSeg).chName = {TSInfo.label};
    segments(iSeg).refName = {TSInfo.refSensor};
    segments(iSeg).samplingRate = [TSInfo.samplingRate];
    segments(iSeg).scale = [TSInfo.resolution];
    segments(iSeg).sampleCount = max(segments(iSeg).samplingRate*segments(iSeg).duration);
end
end

function [eventMarkers] = read_nervus_header_events(h, StaticPackets, Index)
%% Get events  - Andrei Barborica, Dec 2015
% Find sequence of events, that are stored in the section tagged 'Events'
eventsSection = strcmp({StaticPackets.tag}, 'Events');
idxSection = find(eventsSection);
indexIdx = find([Index.sectionIdx] == StaticPackets(idxSection).index);
offset = Index(indexIdx).offset;

ePktLen = 272;    % Event packet length, see EVENTPACKET definition
eMrkLen = 240;    % Event marker length, see EVENTMARKER definition
evtPktGUID = hex2dec({'80', 'F6', '99', 'B7', 'A4', '72', 'D3', '11', '93', 'D3', '00', '50', '04', '00', 'C1', '48'}); % GUID for event packet header
HCEVENT_ANNOTATION        =  '{A5A95612-A7F8-11CF-831A-0800091B5BDA}';
HCEVENT_SEIZURE           =  '{A5A95646-A7F8-11CF-831A-0800091B5BDA}';
HCEVENT_FORMATCHANGE      =  '{08784382-C765-11D3-90CE-00104B6F4F70}';
HCEVENT_PHOTIC            =  '{6FF394DA-D1B8-46DA-B78F-866C67CF02AF}';
HCEVENT_POSTHYPERVENT     =  '{481DFC97-013C-4BC5-A203-871B0375A519}';
HCEVENT_REVIEWPROGRESS    =  '{725798BF-CD1C-4909-B793-6C7864C27AB7}';
HCEVENT_EXAMSTART         =  '{96315D79-5C24-4A65-B334-E31A95088D55}';
HCEVENT_HYPERVENTILATION  =  '{A5A95608-A7F8-11CF-831A-0800091B5BDA}';                            
HCEVENT_IMPEDANCE         =  '{A5A95617-A7F8-11CF-831A-0800091B5BDA}';
DAYSECS = 86400.0;  % From nrvdate.h

fseek(h,offset,'bof');
pktGUID = fread(h,16,'uint8');
pktLen  = fread(h,1,'uint64');
eventMarkers = struct();
i = 0;    % Event counter
while (pktGUID == evtPktGUID)
    i = i + 1;
    % Please refer to EVENTMARKER structure in the Nervus file documentation
    fseek(h,8,'cof'); % Skip eventID, not used
    evtDate           = fread(h,1,'double');
    evtDateFraction   = fread(h,1,'double');
    eventMarkers(i).dateOLE = evtDate;
    eventMarkers(i).dateFraction = evtDateFraction;
    evtPOSIXTime = evtDate*DAYSECS + evtDateFraction - 2209161600;% 2208988800; %8
    eventMarkers(i).dateStr = datestr(evtPOSIXTime/DAYSECS + datenum(1970,1,1),'dd-mmmm-yyyy HH:MM:SS.FFF'); % Save fractions of seconds, as well
    eventMarkers(i).duration  = fread(h,1,'double');
    fseek(h,48,'cof');
    evtUser                       = fread(h,12,'uint16');
    eventMarkers(i).user      = deblank(char(evtUser).');
    evtTextLen                    = fread(h,1,'uint64');
    evtGUID                       = fread(h,16,'uint8');
    eventMarkers(i).GUID      = sprintf('{%.2X%.2X%.2X%.2X-%.2X%.2X-%.2X%.2X-%.2X%.2X-%.2X%.2X%.2X%.2X%.2X%.2X}',evtGUID([4 3 2 1 6 5 8 7 9:16]));
    fseek(h,16,'cof');    % Skip Reserved4 array
    evtLabel                      = fread(h,32,'uint16'); % LABELSIZE = 32;
    evtLabel                      = deblank(char(evtLabel).');    % Not used
    eventMarkers(i).label = evtLabel;
    
    % Only a subset of all event types are dealt with
    switch eventMarkers(i).GUID
        case HCEVENT_SEIZURE
            eventMarkers(i).IDStr = 'Seizure';
            %disp(' Seizure event');
        case HCEVENT_ANNOTATION
            eventMarkers(i).IDStr = 'Annotation';
            fseek(h,32,'cof');    % Skip Reserved5 array
            evtAnnotation = fread(h,evtTextLen,'uint16');
            eventMarkers(i).annotation = deblank(char(evtAnnotation).');
            %disp(sprintf(' Annotation:%s',evtAnnotation));
        case HCEVENT_FORMATCHANGE
            eventMarkers(i).IDStr = 'Format change';
        case HCEVENT_PHOTIC
            eventMarkers(i).IDStr = 'Photic';
        case HCEVENT_POSTHYPERVENT
            eventMarkers(i).IDStr = 'Posthyperventilation';
        case HCEVENT_REVIEWPROGRESS 
            eventMarkers(i).IDStr = 'Review progress';
        case HCEVENT_EXAMSTART
            eventMarkers(i).IDStr = 'Exam start';
        case HCEVENT_HYPERVENTILATION
            eventMarkers(i).IDStr = 'Hyperventilation';
        case HCEVENT_IMPEDANCE
            eventMarkers(i).IDStr = 'Impedance';
        otherwise
            eventMarkers(i).IDStr = 'UNKNOWN';
    end
    
    % Next packet
    offset = offset + pktLen;
    fseek(h,offset,'bof');
    pktGUID = fread(h,16,'uint8');
    pktLen  = fread(h,1,'uint64');
end
end

function [montage] = read_nervus_header_montage(h, StaticPackets, Index)
%% Get montage  - Andrei Barborica, Dec 2015
% Derivation (montage)
mtgIdx  = StaticPackets(find(strcmp({StaticPackets.IDStr},'DERIVATIONGUID'),1)).index;
indexIdx      = find([Index.sectionIdx]==mtgIdx,1);
fseek(h,Index(indexIdx(1)).offset + 40,'bof');    % Beginning of current montage name
mtgName       = deblank(char(fread(h,32,'uint16')).');
fseek(h,640,'cof');                             % Number of traces in the montage
numDerivations = fread(h,1,'uint32');
numDerivations2 = fread(h,1,'uint32');

montage = struct();
for i = 1:numDerivations
    montage(i).derivationName = deblank(char(fread(h,64,'uint16')).');
    montage(i).signalName1    = deblank(char(fread(h,32,'uint16')).');
    montage(i).signalName2    = deblank(char(fread(h,32,'uint16')).');
    fseek(h,264,'cof');         % Skip additional info
end

% Display properties
dispIdx = StaticPackets(find(strcmp({StaticPackets.IDStr},'DISPLAYGUID'),1)).index;
indexIdx  = find([Index.sectionIdx]==dispIdx,1);
fseek(h,Index(indexIdx(1)).offset + 40,'bof');    % Beginning of current montage name
displayName = deblank(char(fread(h,32,'uint16')).');
fseek(h,640,'cof'); % Number of traces in the montage
numTraces = fread(h,1,'uint32');
numTraces2 = fread(h,1,'uint32');

if (numTraces == numDerivations)
    for i = 1:numTraces
        fseek(h,32,'cof');
        montage(i).color = fread(h,1,'uint32'); % Use typecast(uint32(montage(i).color),'uint8') to convert to RGB array
        fseek(h,136-4,'cof');
    end
else
    disp('Could not match montage derivations with display color table');
end
end


function [montages] = read_nervus_header_dynamic_montages(DynamicPackets)
%% Get montages from the dynamic packets  - Jan Brogger - September 2016

montagePackets = DynamicPackets(strcmp({DynamicPackets.IDStr},'DERIVATIONGUID'));
montages = struct();
for numMontage = 1:size(montagePackets,2)
    montagePacket = montagePackets(numMontage);
    
    guidmixed = montagePacket.data(1:16);
    guidnonmixed = [guidmixed(04), guidmixed(03), guidmixed(02), guidmixed(01), ...
        guidmixed(06), guidmixed(05), guidmixed(08), guidmixed(07), ...
        guidmixed(09), guidmixed(10), guidmixed(11), guidmixed(12), ...
        guidmixed(13), guidmixed(15), guidmixed(15), guidmixed(16)];
    montages(numMontage).guid1 = num2str(guidnonmixed,'%02X');
        
    montages(numMontage).packetSize = typecast(montagePacket.data(17:24),'uint64');
    
    guidmixed2 = montagePacket.data(25:40);
    guidnonmixed2 = [guidmixed2(04), guidmixed2(03), guidmixed2(02), guidmixed2(01), ...
        guidmixed2(06), guidmixed2(05), guidmixed2(08), guidmixed2(07), ...
        guidmixed2(09), guidmixed2(10), guidmixed2(11), guidmixed2(12), ...
        guidmixed2(13), guidmixed2(15), guidmixed2(15), guidmixed2(16)];
    montages(numMontage).guid2 = num2str(guidnonmixed2,'%02X');
    
    montages(numMontage).itemName = deblank(cast(montagePacket.data(41:104),'char')');    
    montages(numMontage).elements = typecast(montagePacket.data(745:748),'uint32');
    montages(numMontage).alloc = typecast(montagePacket.data(749:752),'uint32');
    
    offset = 753;
    for i=1:montages(numMontage).elements        
        montages(numMontage).channel(i).name = deblank(cast(montagePacket.data(offset:(offset+127)),'char')');
        offset = offset + 128;
        montages(numMontage).channel(i).active = deblank(cast(montagePacket.data(offset:(offset+63)),'char')');
        offset = offset + 64;
        montages(numMontage).channel(i).reference = deblank(cast(montagePacket.data(offset:(offset+63)),'char')');
        offset = offset + 64;
        montages(numMontage).channel(i).isDerived = montagePacket.data(offset);
        offset = offset + 1;
        montages(numMontage).channel(i).isSpecial = montagePacket.data(offset);
        offset = offset + 1;
        offset = offset + 256;        
        offset = offset + 6;
    end        
end

end
