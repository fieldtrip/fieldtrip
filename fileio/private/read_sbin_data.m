function [trialData] = read_sbin_data(filename, hdr, begtrial, endtrial, chanindx)

% READ_SBIN_DATA reads the data from an EGI segmented simple binary format file
%
% Use as
%   [trialData] = read_sbin_data(filename, hdr, begtrial, endtrial, chanindx)
% with
%   filename       name of the input file
%   hdr            header structure, see FT_READ_HEADER
%   begtrial       first trial to read, mutually exclusive with begsample+endsample
%   endtrial       last trial to read,  mutually exclusive with begsample+endsample
%   chanindx       list with channel indices to read
%
% This function returns a 3-D matrix of size Nchans*Nsamples*Ntrials.
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

fh=fopen_or_error([filename],'r');

version = fread(fh,1,'int32');

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
  version = swapbytes(uint32(version)); %hdr.orig.header_array is already byte-swapped
else
    ft_error('ERROR:  This is not a simple binary file.  Note that NetStation does not successfully directly convert EGIS files to simple binary format.\n');
end;

if bitand(version,1) == 0
    unsegmented = 1;
else
    unsegmented = 0;
end;

precision = bitand(version,6);
NChan=hdr.nChans;
Nevent=hdr.orig.header_array(17);

switch precision
    case 2
        trialLength=2*hdr.nSamples*(hdr.nChans+Nevent)+6;
        dataType='int16';
        dataLength=2;
    case 4
        trialLength=4*hdr.nSamples*(hdr.nChans+Nevent)+6;
        dataType='single';
        dataLength=4;
    case 6
        trialLength=8*hdr.nSamples*(hdr.nChans+Nevent)+6;
        dataType='double';
        dataLength=8;
end

if unsegmented
    status = fseek(fh, 36+Nevent*4, 'bof'); %skip over header
    if status==-1
        ft_error('Failure to skip over header of simple binary file.')
    end;
    status = fseek(fh, ((begtrial-1)*(hdr.nChans+Nevent)*dataLength), 'cof'); %skip previous trials
    if status==-1
        ft_error('Failure to skip over previous trials of simple binary file.')
    end;
    if (hdr.orig.header_array(14))==0 && (hdr.orig.header_array(15) > 1) %epoch-marked simple binary file format
        status = fseek(fh, 30, 'bof'); %skip over header
        if status==-1
            ft_error('Failure to skip over header of simple binary file.')
        end;
        status = fseek(fh, ((begtrial-1)*(hdr.nChans+Nevent)*dataLength), 'cof'); %skip previous trials
        if status==-1
            ft_error('Failure to skip over previous trials of simple binary file.')
        end;
        NSamples    = fread(fh,1,'int32',endian);
        NEvent      = fread(fh,1,'int16',endian);
        for j = 1:Nevent
            EventCodes(j,1:4) = char(fread(fh,[1,4],'char',endian));
        end
        switch precision
            case 2
                [temp,count]    = fread(fh,[NChan+NEvent, NSamples],'int16',endian);
            case 4
                [temp,count]    = fread(fh,[NChan+NEvent, NSamples],'single',endian);
            case 6
                [temp,count]    = fread(fh,[NChan+NEvent, NSamples],'double',endian);
        end
        eventData(:,1:NSamples)  = temp( (NChan+1):(NChan+NEvent), 1:NSamples);
        epocSamples=find(eventData(find(strcmp('epoc',cellstr(EventCodes))),:));
        segmentLengths=diff(epocSamples);
        segmentLengths=[segmentLengths length(eventData)-epocSamples(end)+1];
        trialData=zeros(hdr.nChans,hdr.nSamples,endtrial-begtrial+1);
        for iSegment=1:length(segmentLengths)
            if (iSegment >= begtrial) && (iSegment <= endtrial)
                if iSegment ==1
                    startSample=1;
                else
                    startSample=sum(segmentLengths(1:iSegment-1))+1;
                end;
                trialData(:,1:segmentLengths(iSegment),iSegment-begtrial+1) = temp(1:hdr.nChans,startSample:sum(segmentLengths(1:iSegment))); %data zero-padded to maximum epoch length
            end;
        end;
    else
        nSamples  = endtrial-begtrial+1;    %interpret begtrial and endtrial as sample indices
        [trialData count] = fread(fh, [hdr.nChans+Nevent, nSamples],dataType,endian);
        if count < ((hdr.nChans+Nevent) * nSamples)
            ft_error('Failure to read all samples of simple binary file.')
        end;
    end;
else
    status = fseek(fh, 40+length(hdr.orig.CatLengths)+sum(hdr.orig.CatLengths)+Nevent*4, 'bof'); %skip over header
    if status==-1
        ft_error('Failure to skip over header of simple binary file.')
    end;
    status = fseek(fh, (begtrial-1)*trialLength, 'cof'); %skip over initial segments
    if status==-1
        ft_error('Failure to skip over previous trials of simple binary file.')
    end;
    
    trialData=zeros(hdr.nChans,hdr.nSamples,endtrial-begtrial+1);
    
    for segment=1:(endtrial-begtrial+1)
        status = fseek(fh, 6, 'cof'); %skip over segment info
        if status==-1
            ft_error('Failure to skip over segment info of simple binary file.')
        end;
        
        [temp count] = fread(fh, [(hdr.nChans+Nevent), hdr.nSamples],dataType,endian);
        if count < ((hdr.nChans+Nevent) * hdr.nSamples)
            ft_error('Failure to read all samples of simple binary file.')
        end;
        trialData(:,:,segment) = temp(1:hdr.nChans,:);
    end
end
trialData=trialData(chanindx, :,:);
status = fclose(fh);
if status==-1
    ft_error('Failure to close simple binary file.')
end;
