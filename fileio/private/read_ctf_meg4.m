function [meg] = read_ctf_meg4(fname, hdr, begsample, endsample, chanindx)

% READ_CTF_MEG4 reads specified samples from a CTF continous datafile
% It neglects all trial boundaries as if the data was acquired in
% non-continous mode.
%
% Use as
%   [meg] = read_ctf_meg4(filename, hdr, begsample, endsample, chanindx)
% where
%   filename    name of the datafile, including the .meg4 extension
%   header      with all data information (from read_ctf_meg4)
%   begsample   index of the first sample to read
%   endsample   index of the last sample to read
%   chanindx    index of channels to read (optional, default is all)
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

% this can be used for printing detailled user feedback
fb = false;

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
gain = hdr.gainV(chanindx); % only for selected channels
meg  = raw';            % transpose the raw data
for i=1:size(meg,1)
  meg(i,:) = gain(i)*meg(i,:);
end
