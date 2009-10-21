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

% $Log: read_sbin_events.m,v $
% Revision 1.4  2009/09/04 02:44:49  josdie
% Fixed crash for segmented files due to typo in last revision.
%
% Revision 1.3  2009/04/29 10:55:16  jansch
% incorporated handling of unsegmented files
%
% Revision 1.2  2009/03/11 16:12:34  josdie
% Changed category names to cell variables to better support names with differing lengths.
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.5  2008/12/08 09:36:49  roboos
% added cvs log to the matlab files
%

fid=fopen([filename],'r');
if fid==-1
    error('wrong filename')
end

version		= fread(fid,1,'int32');

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
    error('ERROR:  This is not a simple binary file.  Note that NetStation does not successfully directly convert EGIS files to simple binary format.\n');
end;

if bitand(version,1) == 0
    %error('ERROR:  This is an unsegmented file, which is not supported.\n');
    unsegmented = 1;
else
    unsegmented = 0;
end;

precision = bitand(version,6);
if precision == 0
    error('File precision is not defined.');
end;

%		read header...
year		= fread(fid,1,'int16',endian);
month		= fread(fid,1,'int16',endian);
day			= fread(fid,1,'int16',endian);
hour		= fread(fid,1,'int16',endian);
minute		= fread(fid,1,'int16',endian);
second		= fread(fid,1,'int16',endian);
millisecond = fread(fid,1,'int32',endian);
Samp_Rate	= fread(fid,1,'int16',endian);
NChan		= fread(fid,1,'int16',endian);
Gain 		= fread(fid,1,'int16',endian);
Bits 		= fread(fid,1,'int16',endian);
Range 		= fread(fid,1,'int16',endian);
if unsegmented,
    NumCategors = 0;
    NSegments   = 1;
    NSamples    = fread(fid,1,'int32',endian);
    NEvent      = fread(fid,1,'int16',endian);
    for j = 1:NEvent
        EventCodes(j,1:4) = char(fread(fid,[1,4],'char',endian));
    end
    CateNames   = [];
    CatLengths  = [];
    preBaseline = [];
else
    NumCategors	= fread(fid,1,'int16',endian);
    for j = 1:NumCategors
        CatLengths(j)	= fread(fid,1,'int8',endian);
        for i = 1:CatLengths(j)
            CateNames{j}(i)	= char(fread(fid,1,'char',endian));
        end
    end
    NSegments	= fread(fid,1,'int16',endian);
    NSamples	= fread(fid,1,'int32',endian);			% samples per segment
    NEvent		= fread(fid,1,'int16',endian);			% num events per segment
    EventCodes = [];
    for j = 1:NEvent
        EventCodes(j,1:4)	= char(fread(fid,[1,4],'char',endian));
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

eventData	= zeros(NEvent,NSegments*NSamples);
segHdr      = zeros(NSegments,2);

for j = 1:NSegments*(NSamples/nsmp)
    if unsegmented
        %don't know yet
    else
        %read miniheader per segment
        [segHdr(j,1), count]	= fread(fid, 1,'int16',endian);    %cell
        [segHdr(j,2), count]	= fread(fid, 1,'int32',endian);    %time stamp
    end

    switch precision
        case 2
            [temp,count]	= fread(fid,[NChan+NEvent, nsmp],'int16',endian);
        case 4
            [temp,count]	= fread(fid,[NChan+NEvent, nsmp],'single',endian);
        case 6
            [temp,count]	= fread(fid,[NChan+NEvent, nsmp],'double',endian);
    end
    if (NEvent ~= 0)
        eventData(:,((j-1)*nsmp+1):j*nsmp)	= temp( (NChan+1):(NChan+NEvent), 1:nsmp);
    end
end
fclose(fid);

