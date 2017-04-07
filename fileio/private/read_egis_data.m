function dat = read_egis_data(filename, hdr, begtrial, endtrial, chanindx)

% READ_EGIS_DATA reads the data from an EGI EGIS format file
%
% Use as
%   dat = read_egis_data(filename, hdr, begtrial, endtrial, chanindx);
% where
%   filename       name of the input file
%   hdr            header structure, see FT_READ_HEADER
%   begtrial       first trial to read, mutually exclusive with begsample+endsample
%   endtrial       last trial to read,  mutually exclusive with begsample+endsample
%   chanindx       list with channel indices to read
%
% This function returns a 3-D matrix of size Nchans*Nsamples*Ntrials.
% Note that EGIS session files are defined as always being epoched.
% For session files the trials are organized with the members of each cell grouped
% together.  For average files the "trials" (subjects) are organized with the cells
% also grouped together (e.g., "cell1sub1, cell1sub2, ...).
%_______________________________________________________________________
%
%
% Modified from EGI's EGI Toolbox with permission 2007-06-28 Joseph Dien

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

fh=fopen([filename],'r');
if fh==-1
    error('wrong filename')
end
fclose(fh);

[fhdr,chdr,ename,cnames,fcom,ftext] = read_egis_header(filename);
fh=fopen([filename],'r');
fhdr(1)=fread(fh,1,'int32'); %BytOrd
[str,maxsize,cEndian]=computer;
if fhdr(1)==16909060
    if cEndian == 'B'
        endian = 'ieee-be';
    elseif cEndian == 'L'
        endian = 'ieee-le';
    end;
elseif fhdr(1)==67305985
    if cEndian == 'B'
        endian = 'ieee-le';
    elseif cEndian == 'L'
        endian = 'ieee-be';
    end;
else
    error('This is not an EGIS average file.');
end;

if fhdr(2) == -1
    fileType = 'ave';
elseif fhdr(2) == 3
    fileType = 'ses';
else
    error('This is not an EGIS file.');
end;

dat=zeros(hdr.nChans,hdr.nSamples,endtrial-begtrial+1);

%read to end of file header
status=fseek(fh,(130+(2*fhdr(18))+(4*fhdr(19))),'bof');
status=fseek(fh,2,'cof');
for loop=1:fhdr(18)
  temp=fread(fh,80,'uchar');
  theName=strtok(temp);
  theName=strtok(theName,char(0));
  cnames{loop}=deblank(char(theName))';
  status=fseek(fh,fhdr(24+(loop-1))-80,'cof');
end
fcom=fread(fh,fhdr(20),'uchar');
ftext=fread(fh,fhdr(21),'uchar');
fpad=fread(fh,fhdr(22),'uchar');
status=fseek(fh,-2,'cof');

%read to start of desired data
fseek(fh, ((begtrial-1)*hdr.nChans*hdr.nSamples*2), 'cof');

for segment=1:(endtrial-begtrial+1)
    dat(:,:,segment) = fread(fh, [hdr.nChans, hdr.nSamples],'int16',endian);
end
dat=dat(chanindx, :,:);

if fileType == 'ave'
    dat=dat/fhdr(12); %convert to microvolts
elseif fileType == 'ses'
    dat=dat/5; %convert to microvolts (EGIS sess files created by NetStation use a 5 bins per microvolt scaling)
end;

fclose(fh);
