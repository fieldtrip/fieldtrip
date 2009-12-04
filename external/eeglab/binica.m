% binica() - Run stand-alone binary version of runica() from the
%            Matlab command line. Saves time and memory relative
%            to runica().  If stored in a float file, data are not 
%            read into Matlab, and so may be larger than Matlab
%            can handle owing to memory limitations.
% Usage:
%  >> [wts,sph] = binica( datavar,  'key1', arg1, 'key2', arg2 ...);
% else
%  >> [wts,sph] = binica('datafile', chans, frames, 'key1', arg1, ...);
%
% Inputs:
%   datavar       - (chans,frames) data matrix in the Matlab workspace
%   datafile      - quoted 'filename' of float data file multiplexed by channel
%     channels    -   number of channels in datafile (not needed for datavar)
%     frames      -   number of frames (time points) in datafile (only)
%
% Optional flag,argument pairs:
%   'extended'   - int>=0        [0 default: assume no subgaussian comps]
%                  Search for subgaussian comps: 'extended',1 is recommended
%   'pca'        - int>=0        [0 default: don't reduce data dimension]
%                    NB: 'pca' reduction not recommended unless necessary
%   'sphering'   - 'on'/'off'    first 'sphere' the data {default: 'on'}    
%   'lrate'      - (0<float<<1) starting learning rate {default: 1e-4}
%   'blocksize'  - int>=0        [0 default: heuristic, from data size]
%   'maxsteps'   - int>0         {default: 512}
%   'stop'       - (0<float<<<1) stopping learning rate {default: 1e-7} 
%                    NB: 'stop' <= 1e-7 recommended
%   'weightsin'  - Filename string of inital weight matrix of size
%                  (comps,chans) floats, else a weight matrix variable 
%                  in the current Matlab workspace (copied to a local
%                  .inwts files). You may want to reduce the starting 
%                  'lrate' arg (above) when resuming training, and/or 
%                  to reduce the 'stop' arg (above). By default, binary 
%                  ica begins with the identity matrix after sphering. 
%   'verbose     - 'on'/'off'    {default: 'off'}    
%
% Less frequently used input flags:
%   'posact'     - ('on'/'off') Make maximum value for each comp positive.
%                    NB: 'off' recommended. {default: 'off'} 
%   'annealstep' - (0<float<1)   {default: 0.98}
%   'annealdeg'  - (0<n<360)     {default: 60} 
%   'bias'       - 'on'/'off'    {default: 'on'}    
%   'momentum'   - (0<float<1)   {default: 0 = off]
%
% Outputs:
%   wts          - output weights matrix, size (ncomps,nchans)
%   sph          - output sphere matrix, size (nchans,nchans)
%                  Both files are read from float files left on disk
%   stem         - random integer used in the names of the .sc, .wts, 
%                  .sph, and if requested, .intwts files
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 2000 
%
% See also: runica()

% Calls binary translation of runica() by Sigurd Enghoff

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2000 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% $Log: binica.m,v $
% Revision 1.19  2006/11/07 02:37:18  arno
% add third dim compatibility
%
% Revision 1.18  2006/05/27 03:21:03  toby
% bugs
%
% Revision 1.17  2006/01/29 18:20:02  scott
% help msg -sm
%
% Revision 1.16  2006/01/14 18:04:35  scott
% added 'stem' output -sm
%
% Revision 1.15  2006/01/14 17:37:11  scott
% implemented 'weightsin' flag  -sm
%
% Revision 1.14  2006/01/03 20:05:26  scott
% added 'weightsin' argument (not yet tested), reformed help message -sm
%
% Revision 1.13  2005/03/14 19:43:45  arno
% fixing argument passing
%
% Revision 1.12  2005/03/13 19:03:45  scott
% trying to read flag,arg pairs - FAILS!!!
% Arno/Hilit - PLEASE CHECK LINE 219 and etc. ! -Scott
%
% Revision 1.11  2005/03/13 17:14:40  peter
% comments
%
% Revision 1.10  2004/08/19 19:10:41  arno
% did not change anything
%
% Revision 1.9  2004/02/28 00:00:22  arthur
% edit commandline output text
%
% Revision 1.8  2003/10/25 23:52:05  arno
% remove asking for permission to remove file
%
% Revision 1.7  2003/10/07 18:55:03  arno
% undo changes
%
% Revision 1.6  2003/10/07 18:54:30  arno
% testing
%
% Revision 1.5  2002/11/15 01:44:18  arno
% header for web
%
% Revision 1.4  2002/08/05 18:17:06  arno
% debugging directory finding
%
% Revision 1.3  2002/06/25 02:38:49  scott
% calrified error msgs -sm
%
% Revision 1.2  2002/05/01 18:22:36  arno
% making binica available from everywhere
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 08/07/00 Added warning to update icadefs.m -sm
% 09/08/00 Added tmpint to script, weights and sphere files to avoid
%          conflicts when multiple binica sessions run in the same pwd -sm
% 10/07/00 Fixed bug in reading wts when 'pca' ncomps < nchans -sm
% 07/18/01 replaced var ICA with ICABINARY to try to avoid Matlab6 bug -sm
% 11/06/01 add absolute path of files (lines 157-170 & 198) -ad
% 01-25-02 reformated help & license, added links -ad 
 
function [wts,sph,tmpint] = binica(data,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,var21,var22,var23,var24,var25)

if nargin < 1 | nargin > 25
    more on
    help binica
    more off
    return
end
if size(data,3) > 1, data = reshape(data, size(data,1), size(data,2)*size(data,3) ); end;

icadefs % import ICABINARY and SC
if ~exist('SC')
  fprintf('binica(): You need to update your icadefs file to include ICABINARY and SC.\n')
  return
end;
if exist(SC) ~= 2
  fprintf('binica(): No ica source file ''%s'' is in your Matlab path, check...\n', SC);
  return
else
	SC = which(SC);
	fprintf('binica: using source file ''%s''\n',  SC);
end
if exist(ICABINARY) ~= 2
  fprintf('binica(): ica binary ''%s'' is not in your Matlab path, check\n', ICABINARY);
  return
else
	ICABINARYdir = which(ICABINARY);
	if ~isempty(ICABINARYdir)
		fprintf('binica(): using binary ica file ''%s''\n', ICABINARYdir);
	else
		fprintf('binica(): using binary ica file ''\?/%s''\n', ICABINARY);
	end;
end

%
% select random integer 1-10000 to index the binica data files
% make sure no such script file already exists in the pwd
%
scriptfile = SC;
tmpint = 1000;
while exist(scriptfile)
  tmpints = randperm(10000);
  tmpint = int2str(tmpints(tmpint));
  scriptfile = ['binica' tmpint '.sc'];
end
fprintf('scriptfile = %s\n',scriptfile);

nchans = 0;
tmpdata = [];
if ~isstr(data) % data variable given
  if ~exist('data')
    fprintf('\nbinica(): Variable name data not found.\n');
    return
  end
  nchans = size(data,1);
  nframes = size(data,2);
  tmpdata = ['binica' tmpint '.fdt'];
  floatwrite(data,tmpdata);
  datafile = tmpdata;
  firstarg = 2;

else % data filename given
  if ~exist(data)
    fprintf('\nbinica(): File data not found.\n')
    return
  end
  datafile = data;
  if nargin<3
    fprintf(...
'\nbinica(): Data file name must be followed by chans, frames\n');
    return
  end
  nchans = var2;
  nframes = var3;
  if isstr(nchans) | isstr(nframes)
    fprintf(...
'\nbinica(): chans, frames args must be given after data file name\n');
    return
  end
  firstarg = 4;
end

[flags,args] = read_sc(SC); % read flags and args in master SC file

%
% substitute the flags/args pairs in the .sc file
%
arg = firstarg;
if arg > nargin
   fprintf('binica(): no optional (flag, argument) pairs received.\n');
else
 if (nargin-arg+1)/2 > 1
    fprintf('binica(): processing %d (flag, arg) pairs.\n',(nargin-arg+1)/2);
 else
    fprintf('binica(): processing one (flag, arg) pair.\n');
 end
 while arg <= nargin %%%%%%%%%%%% process flags & args %%%%%%%%%%%%%%%%

  eval(['OPTIONFLAG = var' int2str(arg) ';']); 
  % NB: Found that here Matlab var 'FLAG' is (64,3) why!?!?

  if arg == nargin
    fprintf('\nbinica(): Flag %s needs an argument.\n',OPTIONFLAG)
    return
  end
  eval(['Arg = var' int2str(arg+1) ';']);

  if strcmpi(OPTIONFLAG,'pca')
        ncomps = Arg; % get number of components out for reading wts.
  end

  if strcmpi(OPTIONFLAG,'weightsin')
        wtsin = Arg;
        if exist('wtsin') == 2 % file
           fprintf('   setting %s, %s\n','weightsin',Arg);
        elseif exist('wtsin') == 1 % variable
           if size(wtsin,2) ~= nchans
                fprintf('weightsin variable must be of width %d\n',nchans);
                return
           end
        else
           fprintf('weightsin variable not found.\n');
           return
        end
  end
  arg = arg+2;

  nflags = length(flags);
  for f=1:nflags   % replace SC arg with Arg passed from commandline
    if strcmp(OPTIONFLAG,flags{f})
       args{f} = num2str(Arg);
       fprintf('   setting %s, %s\n',flags{f},args{f});
    end
  end
 end
end

%
% insert name of data files, chans and frames
%
for x=1:length(flags)
  if strcmp(flags{x},'DataFile')
     datafile = [pwd '/' datafile];
     args{x} = datafile;
  elseif strcmp(flags{x},'WeightsOutFile')
     weightsfile = ['binica' tmpint '.wts'];
     weightsfile =  [pwd '/' weightsfile];
     args{x} = weightsfile;
  elseif strcmp(flags{x},'SphereFile')
     spherefile = ['binica' tmpint '.sph'];
     spherefile =  [pwd '/' spherefile];
     args{x} = spherefile;
  elseif strcmp(flags{x},'chans')
     args{x} = int2str(nchans);
  elseif strcmp(flags{x},'frames')
     args{x} = int2str(nframes);
  end
end

%
% write the new .sc file
%
fid = fopen(scriptfile,'w');
for x=1:length(flags)
  fprintf(fid,'%s %s\n',flags{x},args{x});
end
if exist('wtsin') % specify WeightsInfile from 'weightsin' flag, arg
     if exist('wtsin') == 1 % variable
       winfn = [pwd '/binica' tmpint '.inwts'];
       floatwrite(wtsin,winfn);
       fprintf('   saving input weights:\n  ');
       weightsinfile = winfn; % weights in file name
     elseif exist(wtsin) == 2 % file
       weightsinfile = wtsin;
       weightsinfile =  [pwd '/' weightsinfile];
     else
       fprintf('binica(): weightsin file|variable not found.\n');
       return
     end 
    eval(['!ls -l ' weightsinfile]);
    fprintf(fid,'%s %s\n','WeightsInFile',weightsinfile);
end
fclose(fid);
if ~exist(scriptfile)
  fprintf('\nbinica(): ica script file %s not written.\n',...
                                   scriptfile);
  return
end
  
%
% %%%%%%%%%%%%%%%%%%%%%% run binary ica %%%%%%%%%%%%%%%%%%%%%%%%%
%
   fprintf('\nRunning ica from script file %s\n',scriptfile);
   if exist('ncomps')
        fprintf('   Finding %d components.\n',ncomps);
   end
   eval_call = ['!' ICABINARY '<' pwd '/' scriptfile];
   eval(eval_call);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
% read in wts and sph results.
%
if ~exist('ncomps')
  ncomps = nchans;
end

wts = floatread(weightsfile,[ncomps Inf],[],0);
if isempty(wts)
   fprintf('\nbinica(): weight matrix not read.\n')
   return
end
sph = floatread(spherefile,[nchans Inf],[],0);
if isempty(sph)
   fprintf('\nbinica():  sphere matrix not read.\n')
   return
end
fprintf('\nbinary ica files left in pwd:\n');
eval(['!ls -l ' scriptfile ' ' weightsfile ' ' spherefile]);
if exist('wtsin')
   eval(['!ls -l ' weightsinfile]);
end
fprintf('\n');

if isstr(data)
  whos wts sph
else
  whos data wts sph
end

%
% If created by binica(), rm temporary data file
% NOTE: doesn't remove the .sc .wts and .fdt files

if ~isempty(tmpdata)
    eval(['!rm -f ' datafile]);
end

%
%%%%%%%%%%%%%%%%%%% included functions %%%%%%%%%%%%%%%%%%%%%%
%
function sout = rmcomment(s,symb)
     n =1;
     while s(n)~=symb % discard # comments
        n = n+1;
     end
     if n == 1
        sout = [];
     else
       sout = s(1:n-1);
     end
    
function sout = rmspace(s)
       n=1;          % discard leading whitespace
       while n<length(s) & isspace(s(n))
          n = n+1;
       end
       if n<length(s)
          sout = s(n:end);
       else
          sout = [];
       end

function [word,sout] = firstword(s)
       n=1;         
       while n<=length(s) & ~isspace(s(n))
          n = n+1;
       end
       if n>length(s)
         word = [];
         sout = s;
       else
         word = s(1:n-1);
         sout = s(n:end);
       end

function [flags,args] = read_sc(master_sc)
%
% read in the master ica script file SC
%
flags = [];
args  = [];
fid = fopen(master_sc,'r');
if fid < 0
  fprintf('\nbinica(): Master .sc file %s not read!\n',master_sc)
     return
end
%
% read SC file info into flags and args lists
%
s = [];
f = 0; % flag number in file
while isempty(s) | s ~= -1
 s = fgetl(fid);
  if s ~= -1
   if ~isempty(s)
     s = rmcomment(s,'#');
     if ~isempty(s)
       f = f+1;
       s = rmspace(s);
       [w s]=firstword(s);
       if ~isempty(s)
          flags{f} = w;
          s = rmspace(s);
          [w s]=firstword(s);
          args{f} = w;
       end
     end
   end
  end
end 
