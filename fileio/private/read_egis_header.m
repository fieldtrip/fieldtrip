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

% $Log: read_egis_header.m,v $
% Revision 1.2  2009/01/22 19:54:32  josdie
% Cell names were being truncated if they contained a space.  Fixed.
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.3  2008/11/04 08:16:11  roboos
% previous fix did not work, this one should work, thanks to Joe
%
% Revision 1.3  2008/11/03 21:31:00  jdien
% Workaround for fileheader size field overflow for large files
%
% Revision 1.2  2008/09/30 07:47:04  roboos
% replaced all occurences of setstr() with char(), because setstr is deprecated by Matlab
%
% Revision 1.1  2007/07/04 13:22:06  roboos
% initial implementation by Joseph Dien with some corrections by Robert
%

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
