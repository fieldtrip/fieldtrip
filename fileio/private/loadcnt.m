% loadcnt() - Load a Neuroscan continuous signal file.
%
% Usage:
%   >> cnt = loadcnt(file, varargin)
%
% Inputs:
%   filename - name of the file with extension
%
% Optional inputs:
%  't1'         - start at time t1, default 0. Warning, events latency
%                 might be innacurate (this is an open issue).
%  'sample1'    - start at sample1, default 0, overrides t1. Warning,
%                 events latency might be innacurate.
%  'lddur'      - duration of segment to load, default = whole file
%  'ldnsamples' - number of samples to load, default = whole file,
%                 overrides lddur
%  'scale'      - ['on'|'off'] scale data to microvolt (default:'on')
%  'dataformat' - ['int16'|'int32'] default is 'int16' for 16-bit data.
%                 Use 'int32' for 32-bit data.
%  'blockread'  - [integer] by default it is automatically determined
%                 from the file header, though sometimes it finds an
%                 incorect value, so you may want to enter a value manually
%                 here (1 is the most standard value).
%  'memmapfile' - ['memmapfile_name'] use this option if the .cnt file
%                 is too large to read in conventially.  The suffix of
%                 the memmapfile_name must be .fdt.  The memmapfile
%                 functions process files based on their suffix, and an
%                 error will occur if you use a different suffix.
%
% Outputs:
%  cnt          - structure with the continuous data and other informations
%               cnt.header
%               cnt.electloc
%               cnt.data
%               cnt.tag
%
% Authors:   Sean Fitzgibbon, Arnaud Delorme, 2000-
%
% Note: function original name was load_scan41.m
%
% Known limitations:
%  For more see http://www.cnl.salk.edu/~arno/cntload/index.html

% Copyright (C) 2000 Sean Fitzgibbon, <psspf@id.psy.flinders.edu.au>
% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [f,lab,ev2p] = loadcnt(filename,varargin)

if ~isempty(varargin)
  r=struct(varargin{:});
else r = [];
end;

try, r.t1;         catch, r.t1=0; end
try, r.sample1;    catch, r.sample1=[]; end
try, r.lddur;      catch, r.lddur=[]; end
try, r.ldnsamples; catch, r.ldnsamples=[]; end
try, r.scale;      catch, r.scale='on'; end
try, r.blockread;  catch, r.blockread = []; end
try, r.dataformat; catch, r.dataformat = 'auto'; end
try, r.memmapfile; catch, r.memmapfile = ''; end


sizeEvent1 = 8  ; %%% 8  bytes for Event1
sizeEvent2 = 19 ; %%% 19 bytes for Event2
sizeEvent3 = 19 ; %%% 19 bytes for Event3

type='cnt';
if nargin ==1
  scan=0;
end

fid = fopen(filename,'r', 'l');
disp(['Loading file ' filename ' ...'])

h.rev               = fread(fid,12,'char');
h.nextfile          = fread(fid,1,'long');
h.prevfile          = fread(fid,1,'ulong');
h.type              = fread(fid,1,'char');
h.id                = fread(fid,20,'char');
h.oper              = fread(fid,20,'char');
h.doctor            = fread(fid,20,'char');
h.referral          = fread(fid,20,'char');
h.hospital          = fread(fid,20,'char');
h.patient           = fread(fid,20,'char');
h.age               = fread(fid,1,'short');
h.sex               = fread(fid,1,'char');
h.hand              = fread(fid,1,'char');
h.med               = fread(fid,20, 'char');
h.category          = fread(fid,20, 'char');
h.state             = fread(fid,20, 'char');
h.label             = fread(fid,20, 'char');
h.date              = fread(fid,10, 'char');
h.time              = fread(fid,12, 'char');
h.mean_age          = fread(fid,1,'float');
h.stdev             = fread(fid,1,'float');
h.n                 = fread(fid,1,'short');
h.compfile          = fread(fid,38,'char');
h.spectwincomp      = fread(fid,1,'float');
h.meanaccuracy      = fread(fid,1,'float');
h.meanlatency       = fread(fid,1,'float');
h.sortfile          = fread(fid,46,'char');
h.numevents         = fread(fid,1,'int');
h.compoper          = fread(fid,1,'char');
h.avgmode           = fread(fid,1,'char');
h.review            = fread(fid,1,'char');
h.nsweeps           = fread(fid,1,'ushort');
h.compsweeps        = fread(fid,1,'ushort');
h.acceptcnt         = fread(fid,1,'ushort');
h.rejectcnt         = fread(fid,1,'ushort');
h.pnts              = fread(fid,1,'ushort');
h.nchannels         = fread(fid,1,'ushort');
h.avgupdate         = fread(fid,1,'ushort');
h.domain            = fread(fid,1,'char');
h.variance          = fread(fid,1,'char');
h.rate              = fread(fid,1,'ushort'); % A USER CLAIMS THAT SAMPLING RATE CAN BE
h.scale             = fread(fid,1,'double'); % FRACTIONAL IN NEUROSCAN WHICH IS
h.veogcorrect       = fread(fid,1,'char');   % OBVIOUSLY NOT POSSIBLE HERE (BUG 606)
h.heogcorrect       = fread(fid,1,'char');
h.aux1correct       = fread(fid,1,'char');
h.aux2correct       = fread(fid,1,'char');
h.veogtrig          = fread(fid,1,'float');
h.heogtrig          = fread(fid,1,'float');
h.aux1trig          = fread(fid,1,'float');
h.aux2trig          = fread(fid,1,'float');
h.heogchnl          = fread(fid,1,'short');
h.veogchnl          = fread(fid,1,'short');
h.aux1chnl          = fread(fid,1,'short');
h.aux2chnl          = fread(fid,1,'short');
h.veogdir           = fread(fid,1,'char');
h.heogdir           = fread(fid,1,'char');
h.aux1dir           = fread(fid,1,'char');
h.aux2dir           = fread(fid,1,'char');
h.veog_n            = fread(fid,1,'short');
h.heog_n            = fread(fid,1,'short');
h.aux1_n            = fread(fid,1,'short');
h.aux2_n            = fread(fid,1,'short');
h.veogmaxcnt        = fread(fid,1,'short');
h.heogmaxcnt        = fread(fid,1,'short');
h.aux1maxcnt        = fread(fid,1,'short');
h.aux2maxcnt        = fread(fid,1,'short');
h.veogmethod        = fread(fid,1,'char');
h.heogmethod        = fread(fid,1,'char');
h.aux1method        = fread(fid,1,'char');
h.aux2method        = fread(fid,1,'char');
h.ampsensitivity    = fread(fid,1,'float');
h.lowpass           = fread(fid,1,'char');
h.highpass          = fread(fid,1,'char');
h.notch             = fread(fid,1,'char');
h.autoclipadd       = fread(fid,1,'char');
h.baseline          = fread(fid,1,'char');
h.offstart          = fread(fid,1,'float');
h.offstop           = fread(fid,1,'float');
h.reject            = fread(fid,1,'char');
h.rejstart          = fread(fid,1,'float');
h.rejstop           = fread(fid,1,'float');
h.rejmin            = fread(fid,1,'float');
h.rejmax            = fread(fid,1,'float');
h.trigtype          = fread(fid,1,'char');
h.trigval           = fread(fid,1,'float');
h.trigchnl          = fread(fid,1,'char');
h.trigmask          = fread(fid,1,'short');
h.trigisi           = fread(fid,1,'float');
h.trigmin           = fread(fid,1,'float');
h.trigmax           = fread(fid,1,'float');
h.trigdir           = fread(fid,1,'char');
h.autoscale         = fread(fid,1,'char');
h.n2                = fread(fid,1,'short');
h.dir               = fread(fid,1,'char');
h.dispmin           = fread(fid,1,'float');
h.dispmax           = fread(fid,1,'float');
h.xmin              = fread(fid,1,'float');
h.xmax              = fread(fid,1,'float');
h.automin           = fread(fid,1,'float');
h.automax           = fread(fid,1,'float');
h.zmin              = fread(fid,1,'float');
h.zmax              = fread(fid,1,'float');
h.lowcut            = fread(fid,1,'float');
h.highcut           = fread(fid,1,'float');
h.common            = fread(fid,1,'char');
h.savemode          = fread(fid,1,'char');
h.manmode           = fread(fid,1,'char');
h.ref               = fread(fid,10,'char');
h.rectify           = fread(fid,1,'char');
h.displayxmin       = fread(fid,1,'float');
h.displayxmax       = fread(fid,1,'float');
h.phase             = fread(fid,1,'char');
h.screen            = fread(fid,16,'char');
h.calmode           = fread(fid,1,'short');
h.calmethod         = fread(fid,1,'short');
h.calupdate         = fread(fid,1,'short');
h.calbaseline       = fread(fid,1,'short');
h.calsweeps         = fread(fid,1,'short');
h.calattenuator     = fread(fid,1,'float');
h.calpulsevolt      = fread(fid,1,'float');
h.calpulsestart     = fread(fid,1,'float');
h.calpulsestop      = fread(fid,1,'float');
h.calfreq           = fread(fid,1,'float');
h.taskfile          = fread(fid,34,'char');
h.seqfile           = fread(fid,34,'char');
h.spectmethod       = fread(fid,1,'char');
h.spectscaling      = fread(fid,1,'char');
h.spectwindow       = fread(fid,1,'char');
h.spectwinlength    = fread(fid,1,'float');
h.spectorder        = fread(fid,1,'char');
h.notchfilter       = fread(fid,1,'char');
h.headgain          = fread(fid,1,'short');
h.additionalfiles   = fread(fid,1,'int');
h.unused            = fread(fid,5,'char');
h.fspstopmethod     = fread(fid,1,'short');
h.fspstopmode       = fread(fid,1,'short');
h.fspfvalue         = fread(fid,1,'float');
h.fsppoint          = fread(fid,1,'short');
h.fspblocksize      = fread(fid,1,'short');
h.fspp1             = fread(fid,1,'ushort');
h.fspp2             = fread(fid,1,'ushort');
h.fspalpha          = fread(fid,1,'float');
h.fspnoise          = fread(fid,1,'float');
h.fspv1             = fread(fid,1,'short');
h.montage           = fread(fid,40,'char');
h.eventfile         = fread(fid,40,'char');
h.fratio            = fread(fid,1,'float');
h.minor_rev         = fread(fid,1,'char');
h.eegupdate         = fread(fid,1,'short');
h.compressed        = fread(fid,1,'char');
h.xscale            = fread(fid,1,'float');
h.yscale            = fread(fid,1,'float');
h.xsize             = fread(fid,1,'float');
h.ysize             = fread(fid,1,'float');
h.acmode            = fread(fid,1,'char');
h.commonchnl        = fread(fid,1,'uchar');
h.xtics             = fread(fid,1,'char');
h.xrange            = fread(fid,1,'char');
h.ytics             = fread(fid,1,'char');
h.yrange            = fread(fid,1,'char');
h.xscalevalue       = fread(fid,1,'float');
h.xscaleinterval    = fread(fid,1,'float');
h.yscalevalue       = fread(fid,1,'float');
h.yscaleinterval    = fread(fid,1,'float');
h.scaletoolx1       = fread(fid,1,'float');
h.scaletooly1       = fread(fid,1,'float');
h.scaletoolx2       = fread(fid,1,'float');
h.scaletooly2       = fread(fid,1,'float');
h.port              = fread(fid,1,'short');
h.numsamples        = fread(fid,1,'ulong');
h.filterflag        = fread(fid,1,'char');
h.lowcutoff         = fread(fid,1,'float');
h.lowpoles          = fread(fid,1,'short');
h.highcutoff        = fread(fid,1,'float');
h.highpoles         = fread(fid,1,'short');
h.filtertype        = fread(fid,1,'char');
h.filterdomain      = fread(fid,1,'char');
h.snrflag           = fread(fid,1,'char');
h.coherenceflag     = fread(fid,1,'char');
h.continuoustype    = fread(fid,1,'char');
h.eventtablepos     = fread(fid,1,'ulong');
h.continuousseconds = fread(fid,1,'float');
h.channeloffset     = fread(fid,1,'long');
h.autocorrectflag   = fread(fid,1,'char');
h.dcthreshold       = fread(fid,1,'uchar');

for n = 1:h.nchannels
  e(n).lab            = deblank(char(fread(fid,10,'char')'));
  e(n).reference      = fread(fid,1,'char');
  e(n).skip           = fread(fid,1,'char');
  e(n).reject         = fread(fid,1,'char');
  e(n).display        = fread(fid,1,'char');
  e(n).bad            = fread(fid,1,'char');
  e(n).n              = fread(fid,1,'ushort');
  e(n).avg_reference  = fread(fid,1,'char');
  e(n).clipadd        = fread(fid,1,'char');
  e(n).x_coord        = fread(fid,1,'float');
  e(n).y_coord        = fread(fid,1,'float');
  e(n).veog_wt        = fread(fid,1,'float');
  e(n).veog_std       = fread(fid,1,'float');
  e(n).snr            = fread(fid,1,'float');
  e(n).heog_wt        = fread(fid,1,'float');
  e(n).heog_std       = fread(fid,1,'float');
  e(n).baseline       = fread(fid,1,'short');
  e(n).filtered       = fread(fid,1,'char');
  e(n).fsp            = fread(fid,1,'char');
  e(n).aux1_wt        = fread(fid,1,'float');
  e(n).aux1_std       = fread(fid,1,'float');
  e(n).senstivity     = fread(fid,1,'float');
  e(n).gain           = fread(fid,1,'char');
  e(n).hipass         = fread(fid,1,'char');
  e(n).lopass         = fread(fid,1,'char');
  e(n).page           = fread(fid,1,'uchar');
  e(n).size           = fread(fid,1,'uchar');
  e(n).impedance      = fread(fid,1,'uchar');
  e(n).physicalchnl   = fread(fid,1,'uchar');
  e(n).rectify        = fread(fid,1,'char');
  e(n).calib          = fread(fid,1,'float');
end

% finding if 32-bits of 16-bits file
% ----------------------------------
begdata = ftell(fid);
if strcmpi(r.dataformat, 'auto')
  r.dataformat = 'int16';
  if (h.nextfile > 0)
    fseek(fid,h.nextfile+52,'bof');
    is32bit = fread(fid,1,'char');
    if (is32bit == 1)
      r.dataformat = 'int32';
    end;
    fseek(fid,begdata,'bof');
  end;
end;
enddata = h.eventtablepos;   % after data
if strcmpi(r.dataformat, 'int16')
  nums    = (enddata-begdata)/h.nchannels/2;
else nums    = (enddata-begdata)/h.nchannels/4;
end;

% number of sample to read
% ------------------------
if ~isempty(r.sample1)
  r.t1      = r.sample1/h.rate;
else
  r.sample1 = r.t1*h.rate;
end;
if strcmpi(r.dataformat, 'int16')
  startpos = r.t1*h.rate*2*h.nchannels;
else startpos = r.t1*h.rate*4*h.nchannels;
end;
if isempty(r.ldnsamples)
  if ~isempty(r.lddur)
    r.ldnsamples = round(r.lddur*h.rate);
  else r.ldnsamples = nums;
  end;
end;

% FIELDTRIP BUGFIX #1412
% In some cases, the orig.header.numsamples = 0, and the output number of samples is wrong.
% In the previous version of loadcnt.m the orig.header.nums field was used (instead of numsamples), which was changed in r5380 to fix bug #1348.
% This bug (1348) was due to loadcnt.m being updated to the most recent version (from neuroscan), which removed the nums field in favor of using numsamples.
% Below is a workaround for when numsamples is incorrect (bug 1412). The reason is unknown (it looks like a neuroscan data-file specific bug).
% I re-added the nums field to loadcnt.m so that it can be used in ft_read_header.m.
% -roevdmei
h.nums = nums;


% channel offset
% --------------
if ~isempty(r.blockread)
  h.channeloffset = r.blockread;
end;
if h.channeloffset > 1
  fprintf('WARNING: reading data in blocks of %d, if this fails, try using option "''blockread'', 1"\n', ...
    h.channeloffset);
end;

disp('Reading data .....')
if type == 'cnt'
  
  % while (ftell(fid) +1 < h.eventtablepos)
  %d(:,i)=fread(fid,h.nchannels,'int16');
  %end
  fseek(fid, startpos, 0);
  % **** This marks the beginning of the code modified for reading
  % large .cnt files
  
  % Switched to r.memmapfile for continuity.  Check to see if the
  % variable exists.  If it does, then the user has indicated the
  % file is too large to be processed in memory.  If the variable
  % is blank, the file is processed in memory.
  if (~isempty(r.memmapfile))
    % open a file for writing
    foutid = fopen(r.memmapfile, 'w') ;
    
    % This portion of the routine reads in a section of the EEG file
    % and then writes it out to the harddisk.
    samples_left = h.nchannels * r.ldnsamples ;
    
    % the size of the data block to be read is limited to 4M
    % samples.  This equates to 16MB and 32MB of memory for
    % 16 and 32 bit files, respectively.
    data_block = 4000000 ;
    max_rows =  data_block / h.nchannels ;
    
    %warning off ;
    max_written = h.nchannels * uint32(max_rows) ;
    %warning on ;
    
    % This while look tracks the remaining samples.  The
    % data is processed in chunks rather than put into
    % memory whole.
    while (samples_left > 0)
      
      % Check to see if the remaining data is smaller than
      % the general processing block by looking at the
      % remaining number of rows.
      to_read = max_rows ;
      if (data_block > samples_left)
        to_read = samples_left / h.nchannels ;
      end ;
      
      % Read data in a relatively small chunk
      temp_dat = fread(fid, [h.nchannels to_read], r.dataformat) ;
      
      % The data is then scaled using the original routine.
      % In the original routine, the entire data set was scaled
      % after being read in.  For this version, scaling occurs
      % after every chunk is read.
      if strcmpi(r.scale, 'on')
        disp('Scaling data .....')
        %%% scaling to microvolts
        for i=1:h.nchannels
          bas=e(i).baseline;sen=e(i).senstivity;cal=e(i).calib;
          mf=sen*(cal/204.8);
          temp_dat(i,:)=(temp_dat(i,:)-bas).*mf;
        end
      end
      
      % Write out data in float32 form to the file name
      % supplied by the user.
      written = fwrite (foutid, temp_dat, 'float32') ;
      
      if (written ~= max_written)
        samples_left = 0 ;
      else
        samples_left = samples_left - written ;
      end ;
      
    end ;
    
    fclose (foutid) ;
    % Set the dat variable.  This gets used later by other
    % EEGLAB functions.
    dat = r.memmapfile ;
    
    % This variable tracks how the data should be read.
    bReadIntoMemory = false ;
  else
    % The memmapfile variable is empty, read into memory.
    bReadIntoMemory = true ;
  end
  
  % This ends the modifications made to read large files.
  % Everything contained within the following if statement is the
  % original code.
  if (bReadIntoMemory == true)
    if h.channeloffset <= 1
      dat=fread(fid, [h.nchannels Inf], r.dataformat);
      if size(dat,2) < r.ldnsamples
        dat=single(dat);
        r.ldnsamples = size(dat,2);
      else
        dat=single(dat(:,1:r.ldnsamples));
      end;
    else
      h.channeloffset = h.channeloffset/2;
      % reading data in blocks
      dat = zeros( h.nchannels, r.ldnsamples, 'single');
      dat(:, 1:h.channeloffset) = fread(fid, [h.channeloffset h.nchannels], r.dataformat)';
      
      counter = 1;
      while counter*h.channeloffset < r.ldnsamples
        dat(:, counter*h.channeloffset+1:counter*h.channeloffset+h.channeloffset) = ...
          fread(fid, [h.channeloffset h.nchannels], r.dataformat)';
        counter = counter + 1;
      end;
    end ;
    
    % ftell(fid)
    if strcmpi(r.scale, 'on')
      disp('Scaling data .....')
      %%% scaling to microvolts
      for i=1:h.nchannels
        bas=e(i).baseline;sen=e(i).senstivity;cal=e(i).calib;
        mf=sen*(cal/204.8);
        dat(i,:)=(dat(i,:)-bas).*mf;
      end % end for i=1:h.nchannels
    end;  % end if (strcmpi(r.scale, 'on')
  end ;
  
  ET_offset = (double(h.prevfile) * (2^32)) + double(h.eventtablepos);    % prevfile contains high order bits of event table offset, eventtablepos contains the low order bits
  fseek(fid, ET_offset, 'bof');
  
  disp('Reading Event Table...')
  eT.teeg   = fread(fid,1,'uchar');
  eT.size   = fread(fid,1,'ulong');
  eT.offset = fread(fid,1,'ulong');
  
  if eT.teeg==2
    nevents=eT.size/sizeEvent2;
    if nevents > 0
      ev2(nevents).stimtype  = [];
      for i=1:nevents
        ev2(i).stimtype      = fread(fid,1,'ushort');
        ev2(i).keyboard      = fread(fid,1,'char');
        temp                 = fread(fid,1,'uint8');
        ev2(i).keypad_accept = bitand(15,temp);
        ev2(i).accept_ev1    = bitshift(temp,-4);
        ev2(i).offset        = fread(fid,1,'long');
        ev2(i).type          = fread(fid,1,'short');
        ev2(i).code          = fread(fid,1,'short');
        ev2(i).latency       = fread(fid,1,'float');
        ev2(i).epochevent    = fread(fid,1,'char');
        ev2(i).accept        = fread(fid,1,'char');
        ev2(i).accuracy      = fread(fid,1,'char');
      end
    else
      ev2 = [];
    end;
  elseif eT.teeg==3  % type 3 is similar to type 2 except the offset field encodes the global sample frame
    nevents=eT.size/sizeEvent3;
    if nevents > 0
      ev2(nevents).stimtype  = [];
      if r.dataformat == 'int32'
        bytes_per_samp = 4;   % I only have 32 bit data, unable to check whether this is necessary,
      else                      % perhaps there is no type 3 file with 16 bit data
        bytes_per_samp = 2;
      end
      for i=1:nevents
        ev2(i).stimtype      = fread(fid,1,'ushort');
        ev2(i).keyboard      = fread(fid,1,'char');
        temp                 = fread(fid,1,'uint8');
        ev2(i).keypad_accept = bitand(15,temp);
        ev2(i).accept_ev1    = bitshift(temp,-4);
        os                   = fread(fid,1,'ulong');
        ev2(i).offset = os * bytes_per_samp * h.nchannels;
        ev2(i).type          = fread(fid,1,'short');
        ev2(i).code          = fread(fid,1,'short');
        ev2(i).latency       = fread(fid,1,'float');
        ev2(i).epochevent    = fread(fid,1,'char');
        ev2(i).accept        = fread(fid,1,'char');
        ev2(i).accuracy      = fread(fid,1,'char');
      end
    else
      ev2 = [];
    end;
  elseif eT.teeg==1
    nevents=eT.size/sizeEvent1;
    if nevents > 0
      ev2(nevents).stimtype  = [];
      for i=1:nevents
        ev2(i).stimtype      = fread(fid,1,'ushort');
        ev2(i).keyboard      = fread(fid,1,'char');
        
        % modified by Andreas Widmann  2005/05/12  14:15:00
        %ev2(i).keypad_accept = fread(fid,1,'char');
        temp                 = fread(fid,1,'uint8');
        ev2(i).keypad_accept = bitand(15,temp);
        ev2(i).accept_ev1    = bitshift(temp,-4);
        % end modification
        
        ev2(i).offset        = fread(fid,1,'long');
      end;
    else
      ev2 = [];
    end;
  else
    disp('Skipping event table (tag != 1,2,3 ; theoritically impossible)');
    ev2 = [];
  end
  
  
  fseek(fid, -1, 'eof');
  t = fread(fid,'char');
  
  f.header   = h;
  f.electloc = e;
  f.data     = dat;
  f.Teeg     = eT;
  f.event    = ev2;
  f.tag=t;
  % Surgical addition of number of samples
  f.ldnsamples = r.ldnsamples ;
  
  %%%% channels labels
  for i=1:h.nchannels
    plab=sprintf('%c',f.electloc(i).lab);
    if i>1
      lab=str2mat(lab,plab);
    else
      lab=plab;
    end
  end
  
  %%%% to change offest in bytes to points
  if ~isempty(ev2)
    if r.sample1 ~= 0
      fprintf(2,'Warning: events imported with a time shift might be innacurate\n');
    end;
    ev2p=ev2;
    ioff=900+(h.nchannels*75); %% initial offset : header + electordes desc
    if strcmpi(r.dataformat, 'int16')
      for i=1:nevents
        ev2p(i).offset=(ev2p(i).offset-ioff)/(2*h.nchannels) - r.sample1; %% 2 short int end
      end
    else % 32 bits
      for i=1:nevents
        ev2p(i).offset=(ev2p(i).offset-ioff)/(4*h.nchannels) - r.sample1; %% 4 short int end
      end
    end;
    f.event = ev2p;
  end;
  
  frewind(fid);
  fclose(fid);
  
end