function [header_array, CateNames, CatLengths, preBaseline] = read_sbin_header(filename)

% READ_SBIN_HEADER reads the header information from an EGI segmented simple binary format file
%
% Use as
%   [header_array, CateNames, CatLengths, preBaseline] = read_sbin_header(filename)
% with
%   header_array     - differs between versions, read code for details
%   CateNames        - category names
%   CatLengths       - length of category names
%   preBaseline      - number of samples in the baseline prior to the baseline event
% and
%   filename    - the name of the data file
%
% Since there is no unique event code for the segmentation event, and hence the baseline period,
% the first event code in the list will be assumed to be the segmentation event.
% NetStation itself simply ignores possible baseline information when importing simple binary files.
%_______________________________________________________________________
%
%
% Modified from EGI's readEGLY.m with permission 2008-03-31 Joseph Dien
%

% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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
% $Id$

fid=fopen_or_error([filename],'r');

version     = fread(fid,1,'int32');
if isempty(version)
    ft_error('ERROR:  This is not a simple binary file.  Note that NetStation does not successfully directly convert EGIS files to simple binary format.');
end;

%check byteorder
[str,maxsize,cEndian]=computer;
if version < 7
    if cEndian == 'B'
        endian = 'ieee-be';
    elseif cEndian == 'L'
        endian = 'ieee-le';
    end;
elseif (version > 6) && ~bitand(version,6)
    if cEndian == 'B'
        endian = 'ieee-le';
    elseif cEndian == 'L'
        endian = 'ieee-be';
    end;
    version = swapbytes(uint32(version));
else
    ft_error('ERROR:  This is not a simple binary file.  Note that NetStation does not successfully directly convert EGIS files to simple binary format.');
end;

if bitand(version,1) == 0
    unsegmented = 1;
else
    unsegmented = 0;
end;

precision = bitand(version,6);
if precision == 0
    ft_error('File precision is not defined.');
end;

%       read header...
year        = fread(fid,1,'int16',endian);
month       = fread(fid,1,'int16',endian);
day         = fread(fid,1,'int16',endian);
hour        = fread(fid,1,'int16',endian);
minute      = fread(fid,1,'int16',endian);
second      = fread(fid,1,'int16',endian);
millisecond = fread(fid,1,'int32',endian);
Samp_Rate   = fread(fid,1,'int16',endian);
NChan       = fread(fid,1,'int16',endian);
Gain        = fread(fid,1,'int16',endian);
Bits        = fread(fid,1,'int16',endian);
Range       = fread(fid,1,'int16',endian);

epochMarked=0;
if unsegmented,
    NumCategors = 0;
    NSegments   = 1;
    NSamples    = fread(fid,1,'int32',endian);
    NEvent      = fread(fid,1,'int16',endian);
    EventCodes='';
    for j = 1:NEvent
        EventCodes(j,1:4) = char(fread(fid,[1,4],'char',endian));
    end
    CateNames   = [];
    CatLengths  = [];
    preBaseline = 0;
    
    if any(strcmp('epoc',cellstr(EventCodes))) && any(strcmp('tim0',cellstr(EventCodes))) %actually epoch-marked segmented file format
        
        %Note that epoch-marked simple binary file format loses cell information so just the one "cell"
        
        switch precision
            case 2
                [temp,count]    = fread(fid,[NChan+NEvent, NSamples],'int16',endian);
            case 4
                [temp,count]    = fread(fid,[NChan+NEvent, NSamples],'single',endian);
            case 6
                [temp,count]    = fread(fid,[NChan+NEvent, NSamples],'double',endian);
        end
        eventData(:,1:NSamples)  = temp( (NChan+1):(NChan+NEvent), 1:NSamples);
        
        NSegments = length(find(eventData(find(strcmp('epoc',cellstr(EventCodes))),:)));
        epocSamples=find(eventData(find(strcmp('epoc',cellstr(EventCodes))),:));
        segmentLengths=diff(epocSamples);
        segmentLengths=[segmentLengths length(eventData)-epocSamples(end)+1];
        NSamples=max(segmentLengths); % samples per segment, will zero pad out to longest length
        preBaseline=min(find(eventData(find(strcmp('tim0',cellstr(EventCodes))),:)))-min(find(eventData(find(strcmp('epoc',cellstr(EventCodes))),:))); %make assumption all prestimulus durations are the same        
    end;
else
    NumCategors = fread(fid,1,'int16',endian);
    for j = 1:NumCategors
        CatLengths(j)   = fread(fid,1,'int8',endian);
        for i = 1:CatLengths(j)
            CateNames{j}(i) = char(fread(fid,1,'char',endian));
        end
    end
    NSegments   = fread(fid,1,'int16',endian);
    NSamples    = fread(fid,1,'int32',endian);          % samples per segment
    NEvent      = fread(fid,1,'int16',endian);          % num events per segment
    EventCodes = '';
    for j = 1:NEvent
        EventCodes(j,1:4)   = char(fread(fid,[1,4],'char',endian));
    end

    preBaseline=0;
    if NEvent > 0

        for j = 1:NSegments
            [segHdr(j,1), count]    = fread(fid, 1,'int16',endian);    %cell
            [segHdr(j,2), count]    = fread(fid, 1,'int32',endian);    %time stamp
            switch precision
                case 2
                    [temp,count]    = fread(fid,[NChan+NEvent, NSamples],'int16',endian);
                case 4
                    [temp,count]    = fread(fid,[NChan+NEvent, NSamples],'single',endian);
                case 6
                    [temp,count]    = fread(fid,[NChan+NEvent, NSamples],'double',endian);
            end
            eventData(:,((j-1)*NSamples+1):j*NSamples)  = temp( (NChan+1):(NChan+NEvent), 1:NSamples);
        end
    end;
    
    if ~unsegmented
        if NEvent == 0
            theEvent=[];
        elseif NEvent == 1
            %assume this is the segmentation event
            theEvent=find(eventData(1,:)>0);
            theEvent=mod(theEvent(1)-1,NSamples)+1;
        else
            %assume the sample that always has an event is the baseline
            %if more than one, choose the earliest one
            totalEventSamples=mod(find(sum(eventData(:,:))')-1,NSamples)+1;
            baselineCandidates = unique(totalEventSamples);
            counters=hist(totalEventSamples, length(baselineCandidates));
            theEvent=min(baselineCandidates(find(ismember(counters,NSegments))));
        end;
        
        preBaseline=theEvent-1;
        if (preBaseline == -1)
            preBaseline =0;
        end;
        if isempty(preBaseline)
            preBaseline =0;
        end;
    end
end

fclose(fid);

header_array    = double([version year month day hour minute second millisecond Samp_Rate NChan Gain Bits Range NumCategors, NSegments, NSamples, NEvent]);
