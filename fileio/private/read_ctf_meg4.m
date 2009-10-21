function [meg] = read_ctf_meg4(fname, hdr, begsample, endsample, chanindx)

% READ_CTF_MEG4 reads specified samples from a CTF continous datafile
% It neglects all trial boundaries as if the data was acquired in
% non-continous mode.
%
% Use as
%   [meg] = read_ctf_meg4(filename, hdr, begsample, endsample, chanindx)
% where
%   filename	name of the datafile, including the .meg4 extension
%   header      with all data information (from read_ctf_meg4)
%   begsample   index of the first sample to read
%   endsample   index of the last sample to read
%   chanindx	index of channels to read (optional, default is all)
%
% See also READ_CTF_MEG4

% "VSM MedTech Ltd. authorizes the release into public domain under the
% GPL licence of the Matlab source code files "read_ctf_res4.m" and
% "read_ctf_meg4.m" by the authors of said files from the F.C. Donders
% Centre, Nijmegen, The Netherlands."

% Author(s): Jim McKay November 1999
% Last revision: Jim McKay
% Copyright (c) 1999-2000 CTF Systems Inc. All Rights Reserved.
%
% modifications Copyright (C) 2002, Ole Jensen
% modifications Copyright (C) 2003, Robert Oostenveld
%
% $Log: read_ctf_meg4.m,v $
% Revision 1.2  2009/05/07 13:25:16  roboos
% added support for old 64-channel CTF files
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.16  2008/09/30 07:47:04  roboos
% replaced all occurences of setstr() with char(), because setstr is deprecated by Matlab
%
% Revision 1.15  2007/09/11 12:18:23  jansch
% new clean implementation to account for very big datasets > 4GB
%
% Revision 1.14  2005/10/04 15:52:09  roboos
% fixed bug that occured when data is split over more than two files (thanks to flodlan)
%
% Revision 1.13  2005/02/18 13:16:58  roboos
% VSM MedTech Ltd. authorised the release of this code in the public domain
% updated the copyrights, updated the help
%
% Revision 1.12  2004/06/21 19:33:08  roberto
% made 2GB warning dependent on global fb flag
%
% Revision 1.11  2003/07/23 15:02:27  roberto
% added check on valid input for read_ctf_meg4, other changes unknown
%
% Revision 1.10  2003/05/22 09:09:41  roberto
% fixed another bug for >2GB files when selected data in within one trial
%
% Revision 1.7  2003/05/21 13:52:29  roberto
% re-implemented support for >2GB files
% improved checking of input arguments
% fixed bug in chanindx indexing for raw data
%
% Revision 1.6  2003/05/19 15:18:50  roberto
% fixed bugs in memory-efficient reading of continuous data
%
% Revision 1.4  2003/04/17 12:37:41  roberto
% changed error for non-continuous files into warning
%
% Revision 1.3  2003/04/01 06:53:35  roberto
% added support for channel selection
% fixed bug with data allocation over multiple trials
%
% Revision 1.2  2003/03/27 08:30:54  roberto
% fixed bug in reading non-multiplexed trial data
% added error checking
%
% Revision 1.1  2003/03/26 13:34:05  roberto
% new implementation
%

% use global flag for feedback
global fb
if isempty(fb)
  fb = 0;
end

nsmp = hdr.nSamples;
ntrl = hdr.nTrials;
nchn = hdr.nChans;

if begsample<1,              error('cannot read before the start of the data');        end
if endsample>nsmp*ntrl*nchn, error('cannot read beyond the end of the data');          end
if begsample>endsample,      error('cannot read a negative number of samples');        end
if nargin<5,                 chanindx = 1:nchn;                                        end
if isempty(chanindx),        error('no channels were specified for reading CTF data'); end

%open the .meg4 file
fid = fopen(fname,'r','ieee-be');
if fid == -1,
  error('could not open datafile');
end

%check whether it is a known format
CTFformat=char(fread(fid, 8, 'uint8'))';
% This function was written for MEG41RS, but also seems to work for some other formats
if ~strcmp(CTFformat(1,1:7),'MEG41CP') && ~strcmp(CTFformat(1,1:7),'MEG4CPT')
  warning('meg4 format (%s) is not supported for file %s, trying anyway...', CTFformat(1,1:7), fname);
end

%determine size of .meg4 file
fseek(fid, 0, 'eof');
nbytes   = ftell(fid);

%number of trials per 2GB file FIXME assumes constancy across the .meg4 files
ntrlfile = round((nbytes-8)/(4*nchn*nsmp));
%ntrlfile = (nbytes-8)/(4*nchn*nsmp);
openfile = 0;

%determine which trials have to be read
begtrial = ceil(begsample/nsmp);
endtrial = ceil(endsample/nsmp);
trials   = begtrial:endtrial;

%to ensure correct sample handling in the case of multiple trials
sumsmp = [(begtrial-1):(endtrial-1)]*nsmp;
rawbeg = 1;

minchn  = min(chanindx); %to ensure correct channel handling in the case of blockwise reading
maxchn  = max(chanindx);
loopchn = length(trials)==1 || (maxchn-minchn+1)<=nchn; %decide whether to read in blockwise, or the specified samples per channel

raw     = zeros(endsample-begsample+1, length(chanindx)); %allocate memory
for trllop = 1:length(trials)
  trlnr  = trials(trllop);
  filenr = floor(trlnr/(ntrlfile+0.1));

  %ensure that the correct .meg4 file is open
  if filenr~=openfile && filenr>0,
    fclose(fid);
    nextname = sprintf('%s.%d_meg4', fname(1:(end-5)), filenr);
    if fb
      fprintf('data goes beyond 2GB file boundary, continuing with %s\n', nextname);
    end
    fid = fopen(nextname,'r','ieee-be');
    fseek(fid, 0, 'eof');
    openfile = filenr;
  end

  %this is relative to the current datafile
  rawtrl = mod(trlnr-1, ntrlfile) + 1;
  offset = 8 + 4*(rawtrl-1)*nsmp*nchn;

  %begin and endsamples expressed as samples with respect to the current trial
  tmpbeg = max(begsample-sumsmp(trllop), 1);
  tmpend = min(endsample-sumsmp(trllop), nsmp);
  rawend = rawbeg+tmpend-tmpbeg;

  %either read per channel or read entire trialblock and postselect the channels
  if loopchn,
    for chnlop = 1:length(chanindx)
      %this is relative to the current trial
      chanoffset = 4*(chanindx(chnlop)-1)*nsmp;
      sampoffset = 4*(tmpbeg-1);
      fseek(fid, offset+chanoffset+sampoffset, 'bof');
      [tmp, count]               = fread(fid,[tmpend-tmpbeg+1,1],'int32');
      raw(rawbeg:rawend, chnlop) = tmp;
    end
  else
    %this is relative to the current trial
    chanoffset = 4*(minchn-1)*nsmp;
    fseek(fid, offset+chanoffset, 'bof');
    ntmpchn               = maxchn - minchn + 1;
    [tmp, count]          = fread(fid,[nsmp,length(ntmpchn)],'int32');
    selchn                = chanindx - minchn + 1; %relative to the first channel to be read
    raw(rawbeg:rawend, :) = tmp(tmpbeg:tmpend, selchn);
  end
  rawbeg = rawend+1;
end
fclose(fid);

% multiply the dimensionless values with the calibration value
gain = hdr.gainV(chanindx);	% only for selected channels
meg  = raw';			% transpose the raw data
for i=1:size(meg,1)
  meg(i,:) = gain(i)*meg(i,:);
end
