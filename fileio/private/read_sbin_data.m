function [trialData] = read_sbin_data(filename, hdr, begtrial, endtrial, chanindx)

% READ_SBIN_DATA reads the data from an EGI segmented simple binary format file
%
% Use as
%   [trialData] = read_sbin_data(filename, hdr, begtrial, endtrial, chanindx)
% with
%   filename       name of the input file
%   hdr            header structure, see READ_HEADER
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

% $Log: read_sbin_data.m,v $
% Revision 1.2  2009/04/29 10:55:16  jansch
% incorporated handling of unsegmented files
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.4  2008/12/08 09:36:49  roboos
% added cvs log to the matlab files
%

fh=fopen([filename],'r');
if fh==-1
  error('wrong filename')
end

version	= fread(fh,1,'int32');

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
    error('ERROR:  This is not a simple binary file.  Note that NetStation does not successfully directly convert EGIS files to simple binary format.\n');
end;

if bitand(version,1) == 0
    %error('ERROR:  This is an unsegmented file, which is not supported.\n');
    unsegmented = 1;
else
    unsegmented = 0;
end;

precision = bitand(version,6);
Nevents=hdr.orig.header_array(17);

switch precision
    case 2
        trialLength=2*hdr.nSamples*(hdr.nChans+Nevents)+6;
        dataType='int16';
    case 4
        trialLength=4*hdr.nSamples*(hdr.nChans+Nevents)+6;
        dataType='single';
    case 6
        trialLength=8*hdr.nSamples*(hdr.nChans+Nevents)+6;
        dataType='double';
end

if unsegmented
    %interpret begtrial and endtrial as sample indices
    fseek(fh, 36+Nevents*4, 'bof'); %skip over header
    nSamples  = endtrial-begtrial+1;
    trialData = fread(fh, [hdr.nChans+Nevents, nSamples],dataType,endian);
else
    fseek(fh, 40+length(hdr.orig.CatLengths)+sum(hdr.orig.CatLengths)+Nevents*4, 'bof'); %skip over header
    fseek(fh, (begtrial-1)*trialLength, 'cof'); %skip over initial segments

    trialData=zeros(hdr.nChans,hdr.nSamples,endtrial-begtrial+1);

    for segment=1:(endtrial-begtrial+1)
        fseek(fh, 6, 'cof'); %skip over segment info
        temp = fread(fh, [(hdr.nChans+Nevents), hdr.nSamples],dataType,endian);
        trialData(:,:,segment) = temp(1:hdr.nChans,:);
    end
end
trialData=trialData(chanindx, :,:);
fclose(fh);
