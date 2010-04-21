function [fhdr,chdr,ename,cnames,fcom,ftext] = read_egis_header(filename)

% READ_EGIS_HEADER reads the header information from an EGI EGIS format file
%
% Use as
%   [fhdr chdr] = read_egia_header(filename)
% with
%   fhdr        - the file header information
%   chdr        - the cell header information
%   ename       - experiment name
%   cnames      - cell names
%   fcom        - comments
%   ftext       - general text
% and
%   filename    - the name of the data file
%_______________________________________________________________________
%
%
% Modified from EGI's EGI Toolbox with permission 2007-06-28 Joseph Dien

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
% $Id$

fh=fopen([filename],'r');
if fh==-1
  error('wrong filename')
end

%Prep file for reading
fhdr=zeros(1,23);
status=fseek(fh,0,'bof');

%Read in fhdr fields.
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

fhdr(2)=fread(fh,1,'int16',endian); %HdrVer

if (fhdr(2) ~= -1) && (fhdr(2) ~= 3)
  error('This is not an EGIS average file.');
end;

fhdr(3)=fread(fh,1,'uint16',endian); %LHeader
fhdr(4)=fread(fh,1,'uint32',endian); %LData
status=fseek(fh,80,'cof'); %Skip ExptNam
fhdr(5:10)=fread(fh,6,'int16',endian); %RunDate and RunTime

fhdr(11:23)=fread(fh,13,'int16',endian); %Everything else up to LCellHdr

fhdr(24:(24+fhdr(18)-1))=fread(fh,fhdr(18),'uint16',endian); %LCellHdr for each cell

%Read in chdr
chdr=zeros(fhdr(18),5);
status=fseek(fh,(4*fhdr(19)),'cof');
for loop=1:fhdr(18)
  chdr(loop,1)=fread(fh,1,'int16',endian);
  status=fseek(fh,80,'cof');
  chdr(loop, 2:5)=(fread(fh,4,'int16',endian))';
  lspectot=(chdr(loop,2)*(chdr(loop,5)/2));
  lastcol=(6+lspectot-1);
  if lastcol >= 6
    chdr(loop,lastcol)=0;
    chdr(loop,6:lastcol)=(fread(fh,lspectot,'int16',endian))';
  end
end

%Read experiment name
status=fseek(fh,12,'bof');
tempstr=fread(fh,80,'char');
ename=char(tempstr);

%Read cellnames
status=fseek(fh,(130+(2*fhdr(18))+(4*fhdr(19))),'bof');
status=fseek(fh,2,'cof');
for loop=1:fhdr(18)
  temp=fread(fh,80,'char');
  theName=strtok(temp,0);
  cnames{loop}=deblank(char(theName))';
  status=fseek(fh,fhdr(24+(loop-1))-80,'cof');
end

%Read comment
%status=fseek(fh,fhdr(3),'bof');
%status=fseek(fh,-(fhdr(22)+fhdr(21)+fhdr(20)),'cof');
fcom=fread(fh,fhdr(20),'char');

%Read text
%status=fseek(fh,fhdr(3),'bof');
%status=fseek(fh,-(fhdr(22)+fhdr(21)),'cof');
ftext=fread(fh,fhdr(21),'char');

fclose(fh);
