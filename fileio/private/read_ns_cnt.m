% READ_NS_CNT loads header and/or data from a Neuroscan continuous EEG file 
%
% Usage:
%   >> cnt = read_ns_cnt(file, varargin) 
%
% Inputs:
%   filename - name of the file with extension
%
% Optional inputs:
%  't1'         - start at time t1, default 0
%  'sample1'    - start at sample1, default 0, overrides t1
%  'lddur'      - duration of segment to load, default = whole file
%  'ldnsamples' - number of samples to load, default = whole file, 
%                 overrides lddur
%  'blockread'  - size of the blocks to read. Default is 1.
%  'avgref'     - ['yes'|'no'] average reference. Default 'no'.
%  'avrefchan'  - reference channels. Default none.  
%  'format'     - 16 or 32. Default 16.
%
% Outputs:
%  cnt          - structure with the continuous data and other informations
%
% Known limitations: 
%   Initially I couldn't get the continuous data as they would appear 
% using CNTTOASC or CNTTOBIN (www.neuro.com/neuroscan/download.html). 
% We don't have the code for these functions, so we don't really know how
% the raw data is read. After extensive searches, I realized that for
% my continuous CNT files, data was stored in blocks of 40 unsigned short 
% integers for each channel. I couldn't find where this parameter was 
% specified in the header, so I added the option 'blockread' and input 
% the number 40 by hand { cnt = read_ns_cnt('file.cnt', 'blockread', 40) }.
%   By default the size of the block is 1 and this work for most CNT
% files. For more see http://www.cnl.salk.edu/~arno/cntload/index.html    
%
% Authors: Andrew James & Arnaud Delorme 2000-2001 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2000 Andrew James, CERCO, Toulouse, France 
% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: read_ns_cnt.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.10  2008/11/20 12:59:01  roboos
% correct the number of samples in case of 32 bit file
%
% Revision 1.9  2008/11/13 21:20:08  roboos
% read chars as uint8, this solves problem with i18n (2-byte unicode) and recent matlab versions
%
% Revision 1.8  2008/09/30 07:47:04  roboos
% replaced all occurences of setstr() with char(), because setstr is deprecated by Matlab
%
% Revision 1.7  2007/10/31 16:46:13  roboos
% revert to short format as default
%
% Revision 1.6  2007/07/03 16:10:48  roboos
% added a hack for 32 bit neuroscan format (hdr.nsdf=16|32), this should actually be be done using autodetection
%
% Revision 1.5  2005/10/11 08:35:41  roboos
% close file when only header was read (thanks to Durk)
% fixed typo in help
%
% Revision 1.4  2004/06/28 07:36:53  roberto
% changed from DOS to UNIX linefeeds, I am not completely sure whether I made other changes as well
%
% Revision 1.4  2003/04/10 17:56:09  arno
% removing debuging message
%
% Revision 1.3  2003/04/10 17:50:11  arno
% adding error message
%
% Revision 1.2  2002/10/22 23:53:23  arno
% fopen ieee-le for Mac
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%
% original function by a.c.james 2000-2001
% 'blockread' by arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2001

function r=read_ns_cnt(file, varargin)

if nargin < 1
	help read_ns_cnt;
	return;
end;	

% defaults
[datdir,name,ext]=fileparts(file);
if ~isempty(varargin)
	r=struct(varargin{:});
end;

% add defaults
warning off;
try, r.t1; catch, r.t1=0; end
warning on;
try, r.avrefchan; catch, r.avrefchan=[]; end
try, r.blockread; catch, r.blockread=1; end
try, r.avgref;    catch, r.avgref='non'; end

if ~any(file=='.'), file=[file '.cnt']; end

disp(['Loading file ' file ' ...'])

f=fopen(file, 'rb', 'ieee-le');
if f==-1, error([file ' not found']), end

r.filename=file;
r.rev=freadat(f, 0, 20, 'text');
i=find(r.rev==0); try, r.rev(i(1):end)=''; end

r.nchannels=freadat(f, 370, 1, 'ushort');
numsamples=freadat(f, 864, 1, 'long');  % not accurate, see calculation below
samplespos=900 + 75*r.nchannels;
event.tablepos=freadat(f, 886, 1, 'long');
r.nsamples=(event.tablepos - samplespos)/(2*r.nchannels);

%%%%NEW CHANGE TO READ 32 BIT
if     isfield(r, 'format') && r.format == 32
  r.nsamples = r.nsamples / 2;
elseif isfield(r, 'format') && r.format == 16
  r.nsamples = r.nsamples / 1;
else
  r.nsamples = r.nsamples / 1;
end

r.rate=freadat(f, 376, 1, 'ushort');
r.channeloffset=freadat(f, 932, 1, 'long');
r.dt=1/r.rate;
r.scale=freadat(f, 378, 1, 'double');
r.ampsensitivity=freadat(f, 438, 1, 'float');
r.refelectrode=freadat(f, 540, 10, 'text');
if all(r.refelectrode==0), 
   %%disp('No reference electrode set in file, setting to CZ')
   r.refelectrode(1:2)='CZ'; 
end

% reading all parameters
% ----------------------
%a = freadat(f, 0, 1, 'short');
%for i=1:470
%	a = freadat(f, [], 1, 'short');
%	fprintf('offset %3d value %3d\n', i*2, a);
%	%if mod(i, 10) == 0, fprintf('\n'); end;  	
%end;	

% channel parameters
chandat=freadat(f, 900, [75 r.nchannels], 'char');
r.chan.names=char(chandat(1:9,:))';
r.chan.reference=chandat(11,:);
r.chan.gain=chandat(1+63,:);
r.chan.baseline=freadat(f, 900+47, [1 r.nchannels], 'short', 75);
r.chan.sensitivity=freadat(f, 900+59, [1 r.nchannels], 'float', 75);
r.chan.calib=freadat(f, 900+71, [1 r.nchannels], 'float', 75);
r.microvoltscalar=r.chan.sensitivity.*r.chan.calib/204.8;

r.nevent=0;
fseek(f, event.tablepos, 'bof');
r.event.type=fread(f, 1, 'char');
event.size=fread(f, 1, 'long');
%event.offset=fread(f, 1, 'long')

if r.event.type==1
    event.bytes=8;
elseif r.event.type==2
    event.bytes=19;
else
    error('File format error');
end

r.nevent=event.size/event.bytes;
r.event.stimtype=freadat(f, event.tablepos+9, r.nevent, 'short', event.bytes);  % stimtype
r.event.keyboard=freadat(f, event.tablepos+9+2, r.nevent, 'uchar', event.bytes);  % keyboard
r.event.keypadaccept=freadat(f, event.tablepos+9+3, r.nevent, 'uchar', event.bytes);  % keypadaccept
offset=freadat(f, event.tablepos+9+4, r.nevent, 'long', event.bytes);   % offset
r.event.frame=(offset-samplespos)/(r.nchannels*2);  % samplenumber
r.event.time=r.event.frame/r.rate;

try,
   if r.ldheaderonly==1
      fclose(f);
      return
   end
end

try,
   sample1=r.sample1;
   r.t1=r.sample1*r.dt;
catch,
   try 
      startstim=r.startstim;   % startstim = [stimtype occurrence]
      j=find(r.event.stimtype==startstim(1)); 
      if length(startstim)>1
         j=j(startstim(2)); 
      else
         j=j(1); 
      end
      r.t1=r.event.time(j); 
   end
   sample1=round(r.t1/r.dt);   % first sample to read, zero-based
end
startpos=samplespos+sample1*2*r.nchannels;

try, ldnsamples=r.ldnsamples; catch, try, ldnsamples=round(r.lddur/r.dt); catch, ldnsamples=r.nsamples; end, end
try, ldchan=r.ldchan; catch, ldchan=[1:r.nchannels]; end
if ~isempty(ldchan) & ldchan==-1, ldchan=[1:r.nchannels]; end
r.ldchan=ldchan;

% clip events to read window
i=~(sample1<=r.event.frame & r.event.frame<sample1+ldnsamples);
r.nevent=sum(~i);
r.event.stimtype(i)=[];
r.event.keyboard(i)=[];
r.event.keypadaccept(i)=[];
r.event.frame(i)=[];
r.event.time(i)=[];

try, ldraw=r.ldraw; catch, ldraw=0; end;

%%%%NEW CHANGE TO READ 32 BIT
if isfield(r, 'format') && r.format == 32
    df = 'long';
elseif isfield(r, 'format') && r.format == 16
    df = 'short';
else
    % revert to short format as default
    df = 'short';
end

if ~isempty(ldchan)
   if length(ldchan)==r.nchannels
      % all channels

	  if r.blockread == 1	
      	  %dat=freadat(f, startpos, [r.nchannels ldnsamples], 'short');
          dat=freadat(f, startpos, [r.nchannels ldnsamples], df);
 	  else
     	  dat=zeros( length(ldchan), ldnsamples);
      	  %dat(:, 1:r.blockread)=freadat(f, startpos, [r.blockread r.nchannels], 'short')';
          dat(:, 1:r.blockread)=freadat(f, startpos, [r.blockread r.nchannels], df)';

		  counter = 1;	
 		  while counter*r.blockread < ldnsamples
	      	%dat(:, counter*r.blockread+1:counter*r.blockread+r.blockread) = freadat(f, [], [40 r.nchannels], 'short')';
            dat(:, counter*r.blockread+1:counter*r.blockread+r.blockread) = freadat(f, [], [40 r.nchannels], df)';
			counter = counter + 1;
		  end;
	  end;	

      r.dat=zeros( size(dat,2), length(ldchan));
      if ldraw
         %r.dat=int16(dat)';
         r.dat=int32(dat)';
      else
         for j=1:length(ldchan)
            baseline=r.chan.baseline(ldchan(j));
            if baseline==0
               r.dat(:,j)=r.microvoltscalar(ldchan(j))*dat(j,:)';
            else
               r.dat(:,j)=r.microvoltscalar(ldchan(j))*(dat(j,:)-baseline)';
            end
         end
      end
      
      
      if ~isempty(r.avrefchan)
         avref=0;
         for j=1:length(r.avrefchan)
            avref=avref + r.dat(:, r.avrefchan(j));
         end
         avref=1/length(r.avrefchan)*avref;
         
         for j=1:length(r.chan)
            r.dat(:,j)=r.dat(:,j)-avref;
         end
      end
      
   else
      r.dat=zeros(ldnsamples, length(ldchan));
      for j=1:length(ldchan)
         %dat=freadat(f, startpos+(ldchan(j)-1)*2, ldnsamples, 'short', (r.nchannels-1)*2);
         dat=freadat(f, startpos+(ldchan(j)-1)*2, ldnsamples, df, (r.nchannels-1)*2);
         r.dat(:,j)=r.microvoltscalar(ldchan(j))*(dat'-r.chan.baseline(ldchan(j)));
      end
   end
end

fclose(f);
r.dat = r.dat';

% average referencing
% the reference electrode is equal to sum(r,1)/(elec+1)
% -------------------
switch lower(r.avgref)
	case 'yes', 
		r.dat = r.dat-ones(r.nchannels,1)*sum(r.dat,1)/(r.nchannels+1);
end;

disp done
return;

function y=freadat(f, byte, siz, prec, offset)

if nargin<5, 
   skip=0; 
else
   switch prec
   case 'double', s=8;
   case 'float', s=4;
   case 'long', s=4;
   case 'ulong', s=4;
   case 'int16', s=2;
   case 'short', s=2;
   case 'ushort', s=2;
   case 'uchar', s=1;
   case 'char', s=1;
   case 'text', s=1;
   case 'schar', s=1;
   end
   skip=offset-s;
end

if ~isempty(byte)
	fseek(f, byte, 'bof');
end;

if ~strcmp(prec, 'text')
  if ismember(prec,{'char','uchar','schar'})
    prec = 'uint8';
  end
  y=fread(f, siz, prec, skip);
  %y1=fread(f, siz, 'uint8', skip);
  %y2=fread(f, siz, 'uint8', skip);
  %y = y1 + 256*y2;
else
  y=char(fread(f, siz, 'uint8', skip)');
end

