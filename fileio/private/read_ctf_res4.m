function [hdr] = read_ctf_res4(fname)

% READ_CTF_RES4 reads the header in RES4 format from a CTF dataset
%
% Use as
%   [hdr] = read_ctf_res4(filename)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read header information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fname,'r','ieee-be');

% Check if header file exist
if fid == -1
  errMsg = strcat('Could not open header file:',fname);
  error(errMsg);
end

% First 8 bytes contain filetype, check is fileformat is correct.
% This function was written for MEG41RS, but also seems to work for some other formats
CTFformat=char(fread(fid,8,'uint8'))';
if ~strcmp(CTFformat(1,1:7),'MEG41RS') && ~strcmp(CTFformat(1,1:7),'MEG42RS') && ~strcmp(CTFformat(1,1:7),'MEG3RES')
  warning('res4 format (%s) is not supported for file %s, trying anyway...', CTFformat(1,1:7), fname);
end

% Read the initial parameters
appName       = char(fread(fid,256,'uint8'))' ;
dataOrigin    = char(fread(fid,256,'uint8'))' ;
dataDescrip   = char(fread(fid,256,'uint8'))' ;
no_trial_avgd = fread(fid,1,'int16')          ;
data_time     = char(fread(fid,255,'uint8'))';
data_date     = char(fread(fid,255,'uint8'))';

fseek(fid,1288,'bof');
% Read the general recording parameters
no_samples  = fread(fid,1,'int32');
no_channels = fread(fid,1,'int16');
fseek(fid,2,'cof');         % hole of 2 bytes due to improper alignment
sample_rate = fread(fid,1,'double');
epoch       = fread(fid,1,'double');
no_trials   = fread(fid,1,'int16');
fseek(fid,2,'cof');         % hole of 2 bytes due to improper alignment
preTrigpts=fread(fid,1,'int32');

fseek(fid,1360,'bof');
% read in the meg4Filesetup structure
run_name     = char(fread(fid,32,'uint8')');
run_title    = char(fread(fid,256,'uint8')');
instruments  = char(fread(fid,32,'uint8')');
coll_desc    = char(fread(fid,32,'uint8')');
subj_id      = char(fread(fid,32,'uint8')');
operator     = char(fread(fid,32,'uint8')') ;
sensFilename = char(fread(fid,60,'uint8')') ;

% not nececssary to seek, the file pointer is already at the desired location
% fseek(fid,1836,'bof');

% Read in the run description length
rd_len=fread(fid,1,'uint32');
% Go to the run description and read it in
fseek(fid,1844,'bof');
run_desc=char(fread(fid,rd_len,'uint8')');

% read in the filter information
num_filt=fread(fid,1,'uint16');
for fi=0:(num_filt-1),
  %filt_info=fread(fid,18,'uint8');
  filt_freq =fread(fid,1, 'double');
  filt_class=fread(fid,1, 'uint32');
  filt_type =fread(fid,1, 'uint32');
  num_fparm=fread(fid, 1, 'uint16');
  %num_fparm=filt_info(18);
  if num_fparm ~= 0,
    filt_parm=fread(fid,8*num_fparm,'uint8');
  end % if
end % for fi

% Read in the channel names
for i=1:no_channels,
  temp=fread(fid,32,'uint8')';
  temp(find(temp<32 )) = ' ';       % remove non-printable characters
  temp(find(temp>126)) = ' ';       % remove non-printable characters
  endstr = findstr(temp, '-'); temp(endstr:end) = ' ';  % cut off at '-'
  endstr = findstr(temp, ' '); temp(endstr:end) = ' ';  % cut off at ' '
  chan_name(i,:) = char(temp);      % as char array
  chan_label{i}  = deblank(char(temp)); % as cell array
end %for

% pre-allocate some memory space
sensGain = zeros([no_channels,1]);
qGain    = zeros([no_channels,1]);
ioGain   = zeros([no_channels,1]);
sensType = zeros([no_channels,1]);

% Read in the sensor information
fp = ftell(fid);
for chan=1:no_channels,
  fread(fid,1,'uint8');         % Read and ignore 1 byte from enum
  sensType(chan)=fread(fid,1,'uint8');  % Read sensor type
  fread(fid,2,'uint8');         % Read and ignore originalRunNum
  fread(fid,4,'uint8');         % Read and ignore coilShape
  sensGain(chan)=fread(fid,1,'double'); % Read sensor gain in Phi0/Tesla
  qGain(chan)=fread(fid,1,'double');        % Read qxx gain (usually 2^20 for Q20)
  ioGain(chan)=fread(fid,1,'double');       % Read i/o gain of special sensors (usually 1.0)
  ioOffset(chan)=fread(fid,1,'double');
  numCoils(chan)=fread(fid,1,'int16');
  grad_order_no(chan)=fread(fid,1,'int16');
  fread(fid,4,'uint8');

  % read the coil positions and orientations
  for i=1:8
    Chan(chan).coil(i).pos = fread(fid,3,'double')';
    fread(fid,1,'double');
    Chan(chan).coil(i).ori = fread(fid,3,'double')';
    fread(fid,3,'double');
  end

  % read the coil positions and orientations in head coordinates(?)
  for i=1:8
    Chan(chan).coilHC(i).pos = fread(fid,3,'double')';
    fread(fid,1,'double');
    Chan(chan).coilHC(i).ori = fread(fid,3,'double')';
    fread(fid,3,'double');
  end

  % jump to the next sensor info record
  fseek(fid, fp+chan*1328, 'bof');
end % for chan

% close the header file
fclose(fid);

% according to Tom Holroyd, the sensor types are
%
% meg channels are 5, refmag 0, refgrad 1, adcs 18.
% UPPT001 is 11
% UTRG001 is 11
% SCLK01 is 17
% STIM is 11
% SCLK01 is 17
% EEG057 is 9
% ADC06 is 18
% ADC07 is 18
% ADC16 is 18
% V0 is 15

% assign all the variables that should be outputted as header information
hdr.Fs          = sample_rate;
hdr.nChans      = no_channels;
hdr.nSamples    = no_samples;
hdr.nSamplesPre = preTrigpts;
hdr.timeVec     = (1:no_samples)/sample_rate - preTrigpts/sample_rate - 1/sample_rate;
hdr.nTrials     = no_trials;
hdr.gainV       = ioGain./(qGain.*sensGain);
hdr.ioGain      = ioGain;
hdr.qGain       = qGain;
hdr.sensGain    = sensGain;
hdr.sensType    = sensType;

hdr.label       = chan_label(:);
hdr.nameALL     = chan_name;
hdr.Chan        = Chan;
% hdr.rowMEG      = rowMEG;
% hdr.rowEEG      = rowEEG;
% hdr.rowTRIG     = rowTRIG;
% hdr.rowREF      = rowREF;
% hdr.nameEEG     = [];
% hdr.nameMEG     = [];
% hdr.nameEOG     = [];
% hdr.trigV       = [];
% hdr.SwapData    = [];
