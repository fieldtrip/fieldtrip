function [EventCodes, segHdr, eventData] = read_sbin_events(filename)

% READ_SBIN_EVENTS reads the events information from an EGI segmented simple binary format file
%
% Use as
%   [EventCodes, segHdr, eventData] = read_sbin_events(filename)
% with
%   EventCodes      - if NEvent (from header_array) != 0, then array of 4-char event names
%   segHdr          - condition codes and time stamps for each segment
%   eventData       - if NEvent != 0 then event state for each sample, else 'none'
% and
%   filename    - the name of the data file
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

EventCodes = [];
segHdr = [];
eventData = [];

fid=fopen([filename],'r');
if fid==-1
    error('wrong filename')
end

version     = fread(fid,1,'int32');

%check byteorder
[str,maxsize,cEndian]=computer;
if version < 7
    if cEndian == 'B'
        endian = 'ieee-be';
    elseif cEndian == 'L'
        endian = 'ieee-le';
    end
elseif (version > 6) && ~bitand(version,6)
    if cEndian == 'B'
        endian = 'ieee-le';
    elseif cEndian == 'L'
        endian = 'ieee-be';
    end
    version = swapbytes(uint32(version));
else
    error('ERROR:  This is not a simple binary file.  Note that NetStation does not successfully directly convert EGIS files to simple binary format.\n');
end

if bitand(version,1) == 0
    %error('ERROR:  This is an unsegmented file, which is not supported.\n');
    unsegmented = 1;
else
    unsegmented = 0;
end

precision = bitand(version,6);
if precision == 0
    error('File precision is not defined.');
end

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
if unsegmented
    NumCategors = 0;
    NSegments   = 1;
    NSamples    = fread(fid,1,'int32',endian);
    NEvent      = fread(fid,1,'int16',endian);
    for j = 1:NEvent
        EventCodes(j,1:4) = char(fread(fid,[1,4],'char',endian));
    end
    CateNames   = [];
    CatLengths  = [];
    preBaseline = 0;
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
    EventCodes = [];
    for j = 1:NEvent
        EventCodes(j,1:4)   = char(fread(fid,[1,4],'char',endian));
    end
end



%read the actual events
if unsegmented
    %data are multiplexed
    nsmp = 1;
else
    %data are organized in segments
    nsmp = NSamples;
end

segHdr = zeros(NSegments,2);

if unsegmented
    eventData = false(NEvent,NSegments*NSamples);
    switch precision
        case 2
            dataType = 'int16';
            byteSize=2;
        case 4
            dataType = 'single';
            byteSize=4;
        case 6
            dataType = 'double';
            byteSize=8;
    end
    beg_dat = ftell(fid);
    for ii=1:NEvent
        fseek(fid,beg_dat+(NChan+ii-1)*byteSize,'bof');
        eventData(ii,:) = logical(fread(fid, NSamples, ...
            dataType, (NChan+NEvent-1)*byteSize, endian))';
    end
else
    eventData = zeros(NEvent,NSegments*NSamples);
    for j = 1:NSegments
        %read miniheader per segment
        [segHdr(j,1), count]    = fread(fid, 1,'int16',endian);    %cell
        [segHdr(j,2), count]    = fread(fid, 1,'int32',endian);    %time stamp
        
        switch precision
            case 2
                [temp,count]    = fread(fid,[NChan+NEvent, nsmp],'int16',endian);
            case 4
                [temp,count]    = fread(fid,[NChan+NEvent, nsmp],'single',endian);
            case 6
                [temp,count]    = fread(fid,[NChan+NEvent, nsmp],'double',endian);
        end
        if (NEvent ~= 0)
            eventData(:,((j-1)*nsmp+1):j*nsmp)  = temp( (NChan+1):(NChan+NEvent), 1:nsmp);
        end
    end
end
fclose(fid);
