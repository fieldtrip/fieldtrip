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
% $Log: read_ctf_res4.m,v $
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
% Revision 1.17  2008/09/30 07:47:04  roboos
% replaced all occurences of setstr() with char(), because setstr is deprecated by Matlab
%
% Revision 1.16  2008/07/24 08:51:43  roboos
% added the function declaration to the top, which was accidentaly removed by the previous commit
%
% Revision 1.15  2008/07/24 07:22:49  roboos
% replaced fread..char with uint8, solves problem with 16 bit wide characters (thanks to Erick Ortiz)
%
% Revision 1.14  2007/03/07 08:58:46  roboos
% Do not determine the MEG, REF and EEG channels based on the first character of the channel label, removed rowMEG etc from the header. The relevant information is contained in hdr.sensType.
%
% Revision 1.13  2006/03/06 09:41:23  roboos
% changed some |s into ||s
%
% Revision 1.12  2005/07/28 15:12:27  roboos
% fixed bug in hdr.timeVec (thanks to Sanja Kovacevic)
%
% Revision 1.11  2005/05/26 09:58:08  roboos
% removed the construction of grad, that is now done in a separate fieldtrip function (ctf2grad)
%
% Revision 1.10  2005/05/24 07:34:37  roboos
% added sensType to the output header
%
% Revision 1.9  2005/05/23 11:20:17  roboos
% fixed bug for run description which is longer than 256 characters, improved reading of filter information (thanks to Durk Talsma)
%
% Revision 1.8  2005/04/27 06:18:21  roboos
% added support for MEG42RS format, which seems to work (only tested on a single dataset from MIND Institute in Albuquerque, NM, USA)
%
% Revision 1.7  2005/02/18 13:16:58  roboos
% VSM MedTech Ltd. authorised the release of this code in the public domain
% updated the copyrights, updated the help
%
% Revision 1.6  2004/07/02 11:43:32  roboos
% typographic change in comment
%
% Revision 1.5  2004/06/28 07:34:41  roberto
% added some extra low-level output fields to header, added cm units
%
% Revision 1.4  2003/04/01 07:51:53  roberto
% fixed bug in gradiometer channel labels
%
% Revision 1.3  2003/03/28 17:27:27  roberto
% renamed output gradiometer from gradHC to grad
%
% Revision 1.2  2003/03/24 12:34:33  roberto
% minor changes
%
% Revision 1.1  2003/03/13 14:27:09  roberto
% separated read_ctf_ds into res4-header part and data part
% restructured the res4 part, using (incorrect) CTF documentation
%

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
fseek(fid,2,'cof');			% hole of 2 bytes due to improper alignment
sample_rate = fread(fid,1,'double');
epoch       = fread(fid,1,'double');
no_trials   = fread(fid,1,'int16');
fseek(fid,2,'cof');			% hole of 2 bytes due to improper alignment
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
  temp(find(temp<32 )) = ' ';		% remove non-printable characters
  temp(find(temp>126)) = ' ';		% remove non-printable characters
  endstr = findstr(temp, '-'); temp(endstr:end) = ' ';	% cut off at '-'
  endstr = findstr(temp, ' '); temp(endstr:end) = ' ';	% cut off at ' '
  chan_name(i,:) = char(temp);		% as char array
  chan_label{i}  = deblank(char(temp));	% as cell array
end %for

% pre-allocate some memory space
sensGain = zeros([no_channels,1]);
qGain    = zeros([no_channels,1]);
ioGain   = zeros([no_channels,1]);
sensType = zeros([no_channels,1]);

% Read in the sensor information
fp = ftell(fid);
for chan=1:no_channels,
  fread(fid,1,'uint8');			% Read and ignore 1 byte from enum
  sensType(chan)=fread(fid,1,'uint8');	% Read sensor type
  fread(fid,2,'uint8');			% Read and ignore originalRunNum
  fread(fid,4,'uint8');			% Read and ignore coilShape
  sensGain(chan)=fread(fid,1,'double');	% Read sensor gain in Phi0/Tesla
  qGain(chan)=fread(fid,1,'double');		% Read qxx gain (usually 2^20 for Q20)
  ioGain(chan)=fread(fid,1,'double');		% Read i/o gain of special sensors (usually 1.0)
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
