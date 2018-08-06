% runica() - Perform Independent Component Analysis (ICA) decomposition
%            of input data using the logistic infomax ICA algorithm of 
%            Bell & Sejnowski (1995) with the natural gradient feature 
%            of Amari, Cichocki & Yang, or optionally the extended-ICA 
%            algorithm of Lee, Girolami & Sejnowski, with optional PCA 
%            dimension reduction. Annealing based on weight changes is 
%            used to automate the separation process. 
% Usage:
%         >> [weights,sphere] = runica(data); % train using defaults 
%    else
%         >> [weights,sphere,compvars,bias,signs,lrates,activations] ...
%                             = runica(data,'Key1',Value1',...);
% Input:
%    data     = input data (chans,frames*epochs). 
%               Note that if data consists of multiple discontinuous epochs, 
%               each epoch should be separately baseline-zero'd using
%                  >> data = rmbase(data,frames,basevector);
%
% Optional keywords [argument]:
% 'extended'  = [N] perform tanh() "extended-ICA" with sign estimation 
%               N training blocks. If N > 0, automatically estimate the 
%               number of sub-Gaussian sources. If N < 0, fix number of 
%               sub-Gaussian comps to -N [faster than N>0] (default|0 -> off)
% 'pca'       = [N] decompose a principal component     (default -> 0=off)
%               subspace of the data. Value is the number of PCs to retain.
% 'sphering'  = ['on'/'off'] flag sphering of data      (default -> 'on')
% 'weights'   = [W] initial weight matrix               (default -> eye())
%                            (Note: if 'sphering' 'off', default -> spher())
% 'lrate'     = [rate] initial ICA learning rate (<< 1) (default -> heuristic)
% 'block'     = [N] ICA block size (<< datalength)      (default -> heuristic)
% 'anneal'    = annealing constant (0,1] (defaults -> 0.90, or 0.98, extended)
%                         controls speed of convergence
% 'annealdeg' = [N] degrees weight change for annealing (default -> 70)
% 'stop'      = [f] stop training when weight-change < this (default -> 1e-6
%               if less than 33 channel and 1E-7 otherwise)
% 'maxsteps'  = [N] max number of ICA training steps    (default -> 512)
% 'bias'      = ['on'/'off'] perform bias adjustment    (default -> 'on')
% 'momentum'  = [0<f<1] training momentum               (default -> 0)
% 'specgram'  = [srate loHz hiHz frames winframes] decompose a complex time/frequency
%               transform of the data - though not optimally. (Note: winframes must 
%               divide frames) (defaults [srate 0 srate/2 size(data,2) size(data,2)])
% 'posact'    = make all component activations net-positive(default 'off'}
%               Requires time and memory; posact() may be applied separately.
% 'verbose'   = give ascii messages ('on'/'off')        (default -> 'on')
% 'logfile'   = [filename] save all message in a log file in addition to showing them
%               on screen (default -> none)
% 'interput'  = ['on'|'off'] draw interupt figure. Default is off.
%
% Outputs:    [Note: RO means output in reverse order of projected mean variance
%                    unless starting weight matrix passed ('weights' above)]
% weights     = ICA weight matrix (comps,chans)      [RO]
% sphere      = data sphering matrix (chans,chans) = spher(data)
%               Note that unmixing_matrix = weights*sphere {if sphering off -> eye(chans)}
% compvars    = back-projected component variances   [RO]
% bias        = vector of final (ncomps) online bias [RO]    (default = zeros())
% signs       = extended-ICA signs for components    [RO]    (default = ones())
%                   [ -1 = sub-Gaussian; 1 = super-Gaussian]
% lrates      = vector of learning rates used at each training step [RO]
% activations = activation time courses of the output components (ncomps,frames*epochs)
%
% Authors: Scott Makeig with contributions from Tony Bell, Te-Won Lee, 
% Tzyy-Ping Jung, Sigurd Enghoff, Michael Zibulevsky, Delorme Arnaud,
% CNL/The Salk Institute, La Jolla, 1996-

% Uses: posact()

% 'ncomps'    = [N] number of ICA components to compute (default -> chans or 'pca' arg) 
%               using rectangular ICA decomposition. This parameter may return 
%               strange results. This is because the weight matrix is rectangular 
%               instead of being square. Do not use except to try to fix the problem. 

% Reference (please cite):
%
% Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J.,
% "Independent component analysis of electroencephalographic data," 
% In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural 
% Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).
%
% Toolbox Citation:
%
% Makeig, Scott et al. "EEGLAB: ICA Toolbox for Psychophysiological Research". 
% WWW Site, Swartz Center for Computational Neuroscience, Institute of Neural
% Computation, University of San Diego California
% <www.sccn.ucsd.edu/eeglab/>, 2000. [World Wide Web Publication]. 
%
% For more information:
% http://www.sccn.ucsd.edu/eeglab/icafaq.html - FAQ on ICA/EEG
% http://www.sccn.ucsd.edu/eeglab/icabib.html - mss. on ICA & biosignals
% http://www.cnl.salk.edu/~tony/ica.html - math. mss. on ICA

% Copyright (C) 1996 Scott Makeig et al, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit history %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  runica()  - by Scott Makeig with contributions from Tony Bell, Te-Won Lee 
%              Tzyy-Ping Jung, Sigurd Enghoff, Michael Zibulevsky et al.
%                            CNL / Salk Institute 1996-00
%  04-30-96 built from icatest.m and ~jung/.../wtwpwica.m -sm
%  07-28-97 new runica(), adds bias (default on), momentum (default off), 
%           extended-ICA (Lee & Sejnowski, 1997), cumulative angledelta 
%           (until lrate drops), keywords, signcount for speeding extended-ICA
%  10-07-97 put acos() outside verbose loop; verbose 'off' wasn't stopping -sm
%  11-11-97 adjusted help msg -sm
%  11-30-97 return eye(chans) if sphering 'off' or 'none' (undocumented option) -sm
%  02-27-98 use pinv() instead of inv() to rank order comps if ncomps < chans -sm
%  04-28-98 added 'posact' and 'pca' flags  -sm
%  07-16-98 reduced length of randperm() for kurtosis subset calc. -se & sm
%  07-19-98 fixed typo in weights def. above -tl & sm
%  12-21-99 added 'specgram' option suggested by Michael Zibulevsky, UNM -sm
%  12-22-99 fixed rand() sizing inefficiency on suggestion of Mike Spratling, UK -sm
%  01-11-00 fixed rand() sizing bug on suggestion of Jack Foucher, Strasbourg -sm
%  12-18-00 test for existence of Sig Proc Tlbx function 'specgram'; improve
%           'specgram' option arguments -sm
%  01-25-02 reformated help & license -ad 
%  01-25-02 lowered default lrate and block -ad 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [weights,sphere,meanvar,bias,signs,lrates,data,y] = runica(data,varargin) % NB: Now optionally returns activations as variable 'data' -sm 7/05

if nargin < 1
  help runica  
  return
end

[chans frames] = size(data); % determine the data size
urchans = chans;  % remember original data channels 
datalength = frames;
if chans<2 
   fprintf('\nrunica() - data size (%d,%d) too small.\n\n', chans,frames);
   return
end
%
%%%%%%%%%%%%%%%%%%%%%% Declare defaults used below %%%%%%%%%%%%%%%%%%%%%%%%
%
MAX_WEIGHT           = 1e8;       % guess that weights larger than this have blown up
DEFAULT_STOP         = 0.000001;  % stop training if weight changes below this
DEFAULT_ANNEALDEG    = 60;        % when angle change reaches this value,
DEFAULT_ANNEALSTEP   = 0.90;      %     anneal by multiplying lrate by this
DEFAULT_EXTANNEAL    = 0.98;      %     or this if extended-ICA
DEFAULT_MAXSTEPS     = 512;       % ]top training after this many steps 
DEFAULT_MOMENTUM     = 0.0;       % default momentum weight

DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
DEFAULT_RESTART_FAC  = 0.9;       % if weights blowup, restart with lrate
                                  % lower by this factor
MIN_LRATE            = 0.000001;  % if weight blowups make lrate < this, quit
MAX_LRATE            = 0.1;       % guard against uselessly high learning rate
DEFAULT_LRATE        = 0.00065/log(chans); 
                                  % heuristic default - may need adjustment
                                  %   for large or tiny data sets!
% DEFAULT_BLOCK        = floor(sqrt(frames/4));  % heuristic default 
DEFAULT_BLOCK          = ceil(min(5*log(frames),0.3*frames)); % heuristic 
                                  % - may need adjustment!
% Extended-ICA option:
DEFAULT_EXTENDED     = 0;         % default off
DEFAULT_EXTBLOCKS    = 1;         % number of blocks per kurtosis calculation
DEFAULT_NSUB         = 1;         % initial default number of assumed sub-Gaussians
                                  % for extended-ICA
DEFAULT_EXTMOMENTUM  = 0.5;       % momentum term for computing extended-ICA kurtosis
MAX_KURTSIZE         = 6000;      % max points to use in kurtosis calculation
MIN_KURTSIZE         = 2000;      % minimum good kurtosis size (flag warning)
SIGNCOUNT_THRESHOLD  = 25;        % raise extblocks when sign vector unchanged
                                  % after this many steps
SIGNCOUNT_STEP       = 2;         % extblocks increment factor 

DEFAULT_SPHEREFLAG   = 'on';      % use the sphere matrix as the default
                                  %   starting weight matrix
DEFAULT_INTERUPT     = 'off';     % figure interuption
DEFAULT_PCAFLAG      = 'off';     % don't use PCA reduction
DEFAULT_POSACTFLAG   = 'off';     % don't use posact(), to save space -sm 7/05
DEFAULT_VERBOSE      = 1;         % write ascii info to calling screen
DEFAULT_BIASFLAG     = 1;         % default to using bias in the ICA update rule
DEFAULT_RESETRANDOMSEED = true;   % default to reset the random number generator to a 'random state'

%                                 
%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargout < 2, 
    fprintf('runica() - needs at least two output arguments.\n');
    return
end
epochs = 1;							 % do not care how many epochs in data

pcaflag    = DEFAULT_PCAFLAG;
sphering   = DEFAULT_SPHEREFLAG;     % default flags
posactflag = DEFAULT_POSACTFLAG;
verbose    = DEFAULT_VERBOSE;
logfile    = [];

block      = DEFAULT_BLOCK;          % heuristic default - may need adjustment!
lrate      = DEFAULT_LRATE;
annealdeg  = DEFAULT_ANNEALDEG;
annealstep = 0;                      % defaults declared below
nochange   = NaN;
momentum   = DEFAULT_MOMENTUM;
maxsteps   = DEFAULT_MAXSTEPS;

weights    = 0;                      % defaults defined below
ncomps     = chans;
biasflag   = DEFAULT_BIASFLAG;

interupt   = DEFAULT_INTERUPT;
extended   = DEFAULT_EXTENDED;
extblocks  = DEFAULT_EXTBLOCKS;
kurtsize   = MAX_KURTSIZE;
signsbias  = 0.02;                   % bias towards super-Gaussian components
extmomentum= DEFAULT_EXTMOMENTUM;    % exp. average the kurtosis estimates
nsub       = DEFAULT_NSUB;
wts_blowup = 0;                      % flag =1 when weights too large
wts_passed = 0;                      % flag weights passed as argument
reset_randomseed = DEFAULT_RESETRANDOMSEED;

%
%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%
   if (nargin> 1 & rem(nargin,2) == 0)
      fprintf('runica(): Even number of input arguments???')
      return
   end
   for i = 1:2:length(varargin) % for each Keyword
      Keyword = varargin{i};
      Value = varargin{i+1};
      if ~isstr(Keyword)
         fprintf('runica(): keywords must be strings')
         return
      end
      Keyword = lower(Keyword); % convert upper or mixed case to lower

      if strcmp(Keyword,'weights') | strcmp(Keyword,'weight')
         if isstr(Value)
            fprintf(...
      'runica(): weights value must be a weight matrix or sphere')
            return
         else
           weights = Value;
           wts_passed =1;
         end
      elseif strcmp(Keyword,'ncomps')
         if isstr(Value)
            fprintf('runica(): ncomps value must be an integer')
            return
         end
         if ncomps < urchans & ncomps ~= Value
            fprintf('runica(): Use either PCA or ICA dimension reduction');
            return
         end
         fprintf('*****************************************************************************************');
         fprintf('************** WARNING: NCOMPS OPTION OFTEN DOES NOT RETURN ACCURATE RESULTS ************');
         fprintf('************** WARNING: IF YOU FIND THE PROBLEM, PLEASE LET US KNOW          ************');
         fprintf('*****************************************************************************************');
         ncomps = Value;
         if ~ncomps,
            ncomps = chans;
         end
      elseif strcmp(Keyword,'pca') 
         if ncomps < urchans & ncomps ~= Value
            fprintf('runica(): Use either PCA or ICA dimension reduction');
            return
         end
         if isstr(Value)
            fprintf(...
'runica(): pca value should be the number of principal components to retain')
            return
         end
         pcaflag = 'on';
         ncomps = Value;
         if ncomps > chans | ncomps < 1,
            fprintf('runica(): pca value must be in range [1,%d]\n',chans)
            return
         end
         chans = ncomps;
       elseif strcmp(Keyword,'interupt') 
         if ~isstr(Value)
           fprintf('runica(): interupt value must be on or off')
           return
         else 
           Value = lower(Value);
           if ~strcmp(Value,'on') & ~strcmp(Value,'off'),
             fprintf('runica(): interupt value must be on or off')
             return
           end
           interupt = Value;
         end
     elseif strcmp(Keyword,'posact') 
         if ~isstr(Value)
           fprintf('runica(): posact value must be on or off')
           return
         else 
           Value = lower(Value);
           if ~strcmp(Value,'on') & ~strcmp(Value,'off'),
             fprintf('runica(): posact value must be on or off')
             return
           end
           posactflag = Value;
         end
      elseif strcmp(Keyword,'lrate')
         if isstr(Value)
            fprintf('runica(): lrate value must be a number')
            return
         end
         lrate = Value;
         if lrate>MAX_LRATE | lrate <0,
           fprintf('runica(): lrate value is out of bounds'); 
           return
         end
         if ~lrate,
            lrate = DEFAULT_LRATE;
         end
      elseif strcmp(Keyword,'block') | strcmp(Keyword,'blocksize')
         if isstr(Value)
            fprintf('runica(): block size value must be a number')
            return
         end
         block = floor(Value);
         if ~block,
           block = DEFAULT_BLOCK; 
         end
      elseif strcmp(Keyword,'stop') | strcmp(Keyword,'nochange') ...
                    | strcmp(Keyword,'stopping')
         if isstr(Value)
            fprintf('runica(): stop wchange value must be a number')
            return
         end
         nochange = Value;
      elseif strcmp(Keyword,'logfile')
         if ~isstr(Value)
            fprintf('runica(): logfile value must be a string')
            return
         end
         logfile = Value;
      elseif strcmp(Keyword,'maxsteps') | strcmp(Keyword,'steps')
         if isstr(Value)
            fprintf('runica(): maxsteps value must be an integer')
            return
         end
         maxsteps = Value;
         if ~maxsteps,
            maxsteps   = DEFAULT_MAXSTEPS;
         end
         if maxsteps < 0
            fprintf('runica(): maxsteps value (%d) must be a positive integer',maxsteps)
            return
         end
      elseif strcmp(Keyword,'anneal') | strcmp(Keyword,'annealstep')
         if isstr(Value)
            fprintf('runica(): anneal step value (%2.4f) must be a number (0,1)',Value)
            return
         end
         annealstep = Value;
         if annealstep <=0 | annealstep > 1,
            fprintf('runica(): anneal step value (%2.4f) must be (0,1]',annealstep)
            return
         end
      elseif strcmp(Keyword,'annealdeg') | strcmp(Keyword,'degrees')
         if isstr(Value)
            fprintf('runica(): annealdeg value must be a number')
            return
         end
         annealdeg = Value;
         if ~annealdeg,
             annealdeg = DEFAULT_ANNEALDEG;
         elseif annealdeg > 180 | annealdeg < 0
          fprintf('runica(): annealdeg (%3.1f) is out of bounds [0,180]',...
                annealdeg);
          return
                                              
         end
      elseif strcmp(Keyword,'momentum')
         if isstr(Value)
            fprintf('runica(): momentum value must be a number')
            return
         end
         momentum = Value;
         if momentum > 1.0 | momentum < 0
          fprintf('runica(): momentum value is out of bounds [0,1]')
          return
         end
      elseif strcmp(Keyword,'sphering') | strcmp(Keyword,'sphereing') ...
                | strcmp(Keyword,'sphere')
         if ~isstr(Value)
           fprintf('runica(): sphering value must be on, off, or none')
           return
         else 
           Value = lower(Value);
           if ~strcmp(Value,'on') & ~strcmp(Value,'off') & ~strcmp(Value,'none'),
             fprintf('runica(): sphering value must be on or off')
             return
           end
           sphering = Value;
         end
      elseif strcmp(Keyword,'bias')
         if ~isstr(Value)
           fprintf('runica(): bias value must be on or off')
           return
         else 
           Value = lower(Value);
           if strcmp(Value,'on') 
              biasflag = 1;
           elseif strcmp(Value,'off'),
              biasflag = 0;
           else
              fprintf('runica(): bias value must be on or off')
              return
           end
         end
      elseif strcmp(Keyword,'specgram') | strcmp(Keyword,'spec')

         if ~exist('specgram') < 2 % if ~exist or defined workspace variable
           fprintf(...
   'runica(): MATLAB Sig. Proc. Toolbox function "specgram" not found.\n')
           return
         end
         if isstr(Value)
           fprintf('runica(): specgram argument must be a vector')
           return
         end
         srate = Value(1);
         if (srate < 0)
             fprintf('runica(): specgram srate (%4.1f) must be >=0',srate)
             return
           end
         if length(Value)>1
           loHz = Value(2);
           if (loHz < 0 | loHz > srate/2)
             fprintf('runica(): specgram loHz must be >=0 and <= srate/2 (%4.1f)',srate/2)
             return
           end
         else
           loHz = 0; % default
         end
         if length(Value)>2
           hiHz = Value(3);
           if (hiHz < loHz | hiHz > srate/2)
             fprintf('runica(): specgram hiHz must be >=loHz (%4.1f) and <= srate/2 (%4.1f)',loHz,srate/2)
             return
           end
         else
           hiHz = srate/2; % default
         end
         if length(Value)>3
           Hzframes = Value(5);
           if (Hzframes<0 | Hzframes > size(data,2))
             fprintf('runica(): specgram frames must be >=0 and <= data length (%d)',size(data,2))
             return
           end
         else
           Hzframes = size(data,2); % default
         end
         if length(Value)>4
           Hzwinlen = Value(4);
           if rem(Hzframes,Hzwinlen) % if winlen doesn't divide frames
             fprintf('runica(): specgram Hzinc must divide frames (%d)',Hzframes)
             return
           end
         else
           Hzwinlen = Hzframes; % default
         end
         Specgramflag = 1; % set flag to perform specgram()

      elseif strcmp(Keyword,'extended') | strcmp(Keyword,'extend')
         if isstr(Value)
           fprintf('runica(): extended value must be an integer (+/-)')
           return
         else
           extended = 1;      % turn on extended-ICA
           extblocks = fix(Value); % number of blocks per kurt() compute
           if extblocks < 0
                nsub = -1*fix(extblocks);  % fix this many sub-Gauss comps
           elseif ~extblocks,
                extended = 0;             % turn extended-ICA off
           elseif kurtsize>frames,   % length of kurtosis calculation
                kurtsize = frames;
                if kurtsize < MIN_KURTSIZE
                   fprintf(...
   'runica() warning: kurtosis values inexact for << %d points.\n',...
                                                         MIN_KURTSIZE);
                end
           end
         end
      elseif strcmp(Keyword,'verbose') 
         if ~isstr(Value)
            fprintf('runica(): verbose flag value must be on or off')
            return
         elseif strcmp(Value,'on'),
             verbose = 1; 
         elseif strcmp(Value,'off'),
             verbose = 0; 
         else
             fprintf('runica(): verbose flag value must be on or off')
             return
         end
      elseif strcmp(Keyword,'reset_randomseed')
         if ischar(Value)
           if strcmp(Value,'yes')
             reset_randomseed = true;
           elseif strcmp(Value,'no')
             reset_randomseed = false;
           else
             fprintf('runica(): not using the reset_randomseed flag, it should be ''yes'',''no'',0, or 1');
           end
         else
           reset_randomseed = Value;
         end
      else
         fprintf('runica(): unknown flag')
         return
      end
   end

%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize weights, etc. %%%%%%%%%%%%%%%%%%%%%%%%
%
if ~annealstep,
  if ~extended,
    annealstep = DEFAULT_ANNEALSTEP;     % defaults defined above
  else
    annealstep = DEFAULT_EXTANNEAL;       % defaults defined above
  end
end % else use annealstep from commandline

if ~annealdeg, 
    annealdeg  = DEFAULT_ANNEALDEG - momentum*90; % heuristic
    if annealdeg < 0,
        annealdeg = 0;
    end
end
if ncomps >  chans | ncomps < 1
    fprintf('runica(): number of components must be 1 to %d.\n',chans);
    return
end
%
%%%%%%%%%%%%%%%%%%%%% Check keyword values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if frames<chans,
    fprintf('runica(): data length (%d) < data channels (%d)!\n',frames,chans)
    return
elseif block < 2,
    fprintf('runica(): block size %d too small!\n',block)
    return
elseif block > frames, 
    fprintf('runica(): block size exceeds data length!\n');
    return
elseif floor(epochs) ~= epochs,
    fprintf('runica(): data length is not a multiple of the epoch length!\n');
    return
elseif nsub > ncomps
    fprintf('runica(): there can be at most %d sub-Gaussian components!\n',ncomps);
    return
end;

if ~isempty(logfile)
    fid = fopen(logfile, 'w');
    if fid == -1, error('Cannot open logfile for writing'); end;
else
    fid = [];
end;
verb = verbose;

if weights ~= 0,                    % initialize weights
  % starting weights are being passed to runica() from the commandline
    if  chans>ncomps & weights ~=0,
        [r,c]=size(weights);
        if r~=ncomps | c~=chans,
            fprintf('runica(): weight matrix must have %d rows, %d columns.\n', ...
                    chans,ncomps);
            return;
        end
    end
    icaprintf(verb,fid,'Using starting weight matrix named in argument list ...\n');
end;   

% 
% adjust nochange if necessary
%
if isnan(nochange) 
    if ncomps > 32
        nochange = 1E-7;
        nochangeupdated = 1; % for fprinting purposes
    else
        nochangeupdated = 1; % for fprinting purposes
        nochange = DEFAULT_STOP;
    end;
else 
    nochangeupdated = 0;
end;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process the data %%%%%%%%%%%%%%%%%%%%%%%%%%
%
icaprintf(verb,fid,'\nInput data size [%d,%d] = %d channels, %d frames\n', ...
          chans,frames,chans,frames);

if strcmp(pcaflag,'on')
    icaprintf(verb,fid,'After PCA dimension reduction,\n  finding ');
else
    icaprintf(verb,fid,'Finding ');
end
if ~extended
    icaprintf(verb,fid,'%d ICA components using logistic ICA.\n',ncomps);
else % if extended
icaprintf(verb,fid,'%d ICA components using extended ICA.\n',ncomps);
if extblocks > 0
    icaprintf(verb,fid,'Kurtosis will be calculated initially every %d blocks using %d data points.\n',...
              extblocks, kurtsize);
else
    icaprintf(verb,fid,'Kurtosis will not be calculated. Exactly %d sub-Gaussian components assumed.\n',nsub);
end
end
icaprintf(verb,fid,'Decomposing %d frames per ICA weight ((%d)^2 = %d weights, %d frames)\n',...
          floor(frames/ncomps.^2),ncomps.^2,frames);
icaprintf(verb,fid,'Initial learning rate will be %g, block size %d.\n',...
          lrate,block);
if momentum>0,
    icaprintf(verb,fid,'Momentum will be %g.\n',momentum);
end
icaprintf(verb,fid,'Learning rate will be multiplied by %g whenever angledelta >= %g deg.\n', ...
          annealstep,annealdeg);

if nochangeupdated 
    icaprintf(verb,fid,'More than 32 channels: default stopping weight change 1E-7\n');
end;
icaprintf(verb,fid,'Training will end when wchange < %g or after %d steps.\n', nochange,maxsteps);
if biasflag,
    icaprintf(verb,fid,'Online bias adjustment will be used.\n');
else
    icaprintf(verb,fid,'Online bias adjustment will not be used.\n');
end

%
%%%%%%%%%%%%%%%%% Remove overall row means of data %%%%%%%%%%%%%%%%%%%%%%%
%
icaprintf(verb,fid,'Removing mean of each channel ...\n');

%BLGBLGBLG replaced
% rowmeans = mean(data');
% data = data - rowmeans'*ones(1,frames);      % subtract row means
%BLGBLGBLG replacement starts
rowmeans = mean(data,2)'; %BLG
% data = data - rowmeans'*ones(1,frames);      % subtract row means
for iii=1:size(data,1) %avoids memory errors BLG
    data(iii,:)=data(iii,:)-rowmeans(iii);
end
%BLGBLGBLG replacement ends

icaprintf(verb,fid,'Final training data range: %g to %g\n', min(min(data)),max(max(data)));

%
%%%%%%%%%%%%%%%%%%% Perform PCA reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if strcmp(pcaflag,'on')
    icaprintf(verb,fid,'Reducing the data to %d principal dimensions...\n',ncomps);
    
    %BLGBLGBLG replaced
    %[eigenvectors,eigenvalues,data] = pcsquash(data,ncomps);
    % make data its projection onto the ncomps-dim principal subspace
    %BLGBLGBLG replacement starts
    %[eigenvectors,eigenvalues,data] = pcsquash(data,ncomps);
    % no need to re-subtract row-means, it was done a few lines above!
    PCdat2 = data';                    % transpose data
    [PCn,PCp]=size(PCdat2);                  % now p chans,n time points
    PCdat2=PCdat2/PCn;
    PCout=data*PCdat2;
    clear PCdat2;
    
    [PCV,PCD] = eig(PCout);                  % get eigenvectors/eigenvalues
    [PCeigenval,PCindex] = sort(diag(PCD));
    PCindex=rot90(rot90(PCindex));
    PCEigenValues=rot90(rot90(PCeigenval))';
    PCEigenVectors=PCV(:,PCindex);
    %PCCompressed = PCEigenVectors(:,1:ncomps)'*data;
    data = PCEigenVectors(:,1:ncomps)'*data;
    
    eigenvectors=PCEigenVectors;
    eigenvalues=PCEigenValues; %#ok<NASGU>
    
    clear PCn PCp PCout PCV PCD PCeigenval PCindex PCEigenValues PCEigenVectors
    %BLGBLGBLG replacement ends
    
end

%
%%%%%%%%%%%%%%%%%%% Perform specgram transformation %%%%%%%%%%%%%%%%%%%%%%%
%
if exist('Specgramflag') == 1
  % [P F T] = SPECGRAM(A,NFFT,Fs,WINDOW,NOVERLAP) % MATLAB Sig Proc Toolbox
  % Hzwinlen =  fix(srate/Hzinc); % CHANGED FROM THIS 12/18/00 -sm

  Hzfftlen = 2^(ceil(log(Hzwinlen)/log(2)));   % make FFT length next higher 2^k
  Hzoverlap = 0; % use sequential windows
  %
  % Get freqs and times from 1st channel analysis
  %
  [tmp,freqs,tms] = specgram(data(1,:),Hzfftlen,srate,Hzwinlen,Hzoverlap);

  fs = find(freqs>=loHz & freqs <= hiHz);
  icaprintf(verb,fid,'runica(): specified frequency range too narrow, exiting!\n');
    
  specdata = reshape(tmp(fs,:),1,length(fs)*size(tmp,2));
  specdata = [real(specdata) imag(specdata)];
     % fprintf('   size(fs) = %d,%d\n',size(fs,1),size(fs,2));
     % fprintf('   size(tmp) = %d,%d\n',size(tmp,1),size(tmp,2));
  %
  % Loop through remaining channels
  %
  for ch=2:chans
      [tmp] = specgram(data(ch,:),Hzwinlen,srate,Hzwinlen,Hzoverlap);
      tmp = reshape((tmp(fs,:)),1,length(fs)*size(tmp,2));
      specdata = [specdata;[real(tmp) imag(tmp)]]; % channels are rows
  end
  %
  % Print specgram confirmation and details
  %
  icaprintf(verb,fid,'Converted data to %d channels by %d=2*%dx%d points spectrogram data.\n',...
            chans,2*length(fs)*length(tms),length(fs),length(tms));
  if length(fs) > 1
      icaprintf(verb,fid,'   Low Hz %g, high Hz %g, Hz incr %g, window length %d\n',freqs(fs(1)),freqs(fs(end)),freqs(fs(2))-freqs(fs(1)),Hzwinlen);
  else
      icaprintf(verb,fid,'   Low Hz %g, high Hz %g, window length %d\n',freqs(fs(1)),freqs(fs(end)),Hzwinlen);
  end
  %
  % Replace data with specdata
  %
  data = specdata;
  datalength=size(data,2);
end
%
%%%%%%%%%%%%%%%%%%% Perform sphering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if strcmp(sphering,'on'), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  icaprintf(verb,fid,'Computing the sphering matrix...\n');
  sphere = 2.0*inv(sqrtm(double(cov(data')))); % find the "sphering" matrix = spher()
  if ~weights,
      icaprintf(verb,fid,'Starting weights are the identity matrix ...\n');
      weights = eye(ncomps,chans); % begin with the identity matrix
  else % weights given on commandline
      icaprintf(verb,fid,'Using starting weights named on commandline ...\n');
  end
  icaprintf(verb,fid,'Sphering the data ...\n');
  data = sphere*data; % decorrelate the electrode signals by 'sphereing' them

elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~weights % is starting weights not given
      icaprintf(verb,fid,'Using the sphering matrix as the starting weight matrix ...\n');
      icaprintf(verb,fid,'Returning the identity matrix in variable "sphere" ...\n');
      sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
      weights = eye(ncomps,chans)*sphere;  % begin with the identity matrix
      sphere = eye(chans);                 % return the identity matrix
  else % weights ~= 0
      icaprintf(verb,fid,'Using starting weights from commandline ...\n');
      icaprintf(verb,fid,'Returning the identity matrix in variable "sphere" ...\n');
      sphere = eye(chans);                 % return the identity matrix
  end
elseif strcmp(sphering,'none')
  sphere = eye(chans,chans);% return the identity matrix
  if ~weights
      icaprintf(verb,fid,'Starting weights are the identity matrix ...\n');
      icaprintf(verb,fid,'Returning the identity matrix in variable "sphere" ...\n');
      weights = eye(ncomps,chans); % begin with the identity matrix
  else % weights ~= 0
      icaprintf(verb,fid,'Using starting weights named on commandline ...\n');
      icaprintf(verb,fid,'Returning the identity matrix in variable "sphere" ...\n');
  end
  icaprintf(verb,fid,'Returned variable "sphere" will be the identity matrix.\n');
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
%
lastt=fix((datalength/block-1)*block+1);
BI=block*eye(ncomps,ncomps);
delta=zeros(1,chans*ncomps);
changes = [];
degconst = 180./pi;
startweights = weights;
prevweights = startweights;
oldweights = startweights;
prevwtchange = zeros(chans,ncomps);
oldwtchange = zeros(chans,ncomps);
lrates = zeros(1,maxsteps);
onesrow = ones(1,block);
bias = zeros(ncomps,1);
signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1
for k=1:nsub
    signs(k) = -1;
end
if extended & extblocks < 0,
    icaprintf(verb,fid,'Fixed extended-ICA sign assignments:  ');
    for k=1:ncomps
        icaprintf(verb,fid,'%d ',signs(k));
    end; icaprintf(verb,fid,'\n');
end
signs = diag(signs); % make a diagonal matrix
oldsigns = zeros(size(signs));
signcount = 0;              % counter for same-signs
signcounts = [];
urextblocks = extblocks;    % original value, for resets
old_kk = zeros(1,ncomps);   % for kurtosis momemtum

%
%%%%%%%% ICA training loop using the logistic sigmoid %%%%%%%%%%%%%%%%%%%
%
icaprintf(verb,fid,'Beginning ICA training ...');
if extended,
    icaprintf(verb,fid,' first training step may be slow ...\n');
else
    icaprintf(verb,fid,'\n');
end
step=0;
laststep=0;
blockno = 1;  % running block counter for kurtosis interrupts

if reset_randomseed
    rand('state',sum(100*clock));  % set the random number generator state to
end                                % a position dependent on the system clock

% interupt figure
% --------------- 
if strcmpi(interupt, 'on')
    fig = figure('visible', 'off');
    supergui( fig, {1 1}, [], {'style' 'text' 'string' 'Press button to interrupt runica()' }, ...
              {'style' 'pushbutton' 'string' 'Interupt' 'callback' 'setappdata(gcf, ''run'', 0);' } );
    set(fig, 'visible', 'on');
    setappdata(gcf, 'run', 1);
    
    if strcmpi(interupt, 'on')
        drawnow;
    end;
end;


%% Compute ICA Weights
if biasflag & extended
    while step < maxsteps, %%% ICA step = pass through all the data %%%%%%%%%
        timeperm=randperm(datalength); % shuffle data order at each step

        for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
            if strcmpi(interupt, 'on')
                drawnow;
                flag = getappdata(fig, 'run');
                if ~flag,
                    if ~isempty(fid), fclose(fid); end;
                    close; error('USER ABORT');
                end;
            end;
            
            %% promote data block (only) to double to keep u and weights double
            u=weights*double(data(:,timeperm(t:t+block-1))) + bias*onesrow;

            y=tanh(u);                                                       
            weights = weights + lrate*(BI-signs*y*u'-u*u')*weights;
            bias = bias + lrate*sum((-2*y)')';  % for tanh() nonlin.
            
            if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                weights = weights + momentum*prevwtchange;
                prevwtchange = weights-prevweights;
                prevweights = weights;
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if max(max(abs(weights))) > MAX_WEIGHT
                wts_blowup = 1;
                change = nochange;
            end
            if ~wts_blowup
                %
                %%%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
                %while step < maxsteps
                if extblocks > 0 & rem(blockno,extblocks) == 0,
                    % recompute signs vector using kurtosis
                    if kurtsize < frames % 12-22-99 rand() size suggestion by M. Spratling
                        rp = fix(rand(1,kurtsize)*datalength);  % pick random subset
                        % Accout for the possibility of a 0 generation by rand
                        ou = find(rp == 0);
                        while ~isempty(ou) % 1-11-00 suggestion by J. Foucher
                            rp(ou) = fix(rand(1,length(ou))*datalength);
                            ou = find(rp == 0);
                        end
                        partact=weights*double(data(:,rp(1:kurtsize)));
                    else                                        % for small data sets,
                        partact=weights*double(data);           % use whole data
                    end
                    m2=mean(partact'.^2).^2;
                    m4= mean(partact'.^4);
                    kk= (m4./m2)-3.0;                           % kurtosis estimates
                    if extmomentum
                        kk = extmomentum*old_kk + (1.0-extmomentum)*kk; % use momentum
                        old_kk = kk;
                    end
                    signs=diag(sign(kk+signsbias));             % pick component signs
                    if signs == oldsigns,
                        signcount = signcount+1;
                    else
                        signcount = 0;
                    end
                    oldsigns = signs;
                    signcounts = [signcounts signcount];
                    if signcount >= SIGNCOUNT_THRESHOLD,
                        extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
                        signcount = 0;                             % less frequent if sign
                    end                                         % is not changing
                end % extblocks > 0 & . . .
            end % if extended & ~wts_blowup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            blockno = blockno + 1;
            if wts_blowup
                break
            end
        end % for t=1:block:lastt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~wts_blowup
            oldwtchange = weights-oldweights;
            step=step+1;
            %
            %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
            %
            lrates(1,step) = lrate;
            angledelta=0.;
            delta=reshape(oldwtchange,1,chans*ncomps);
            change=delta*delta';
        end
        %
        %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
        %
        if wts_blowup | isnan(change)|isinf(change),  % if weights blow up,
            icaprintf(verb,fid,'');
            step = 0;                          % start again
            change = nochange;
            wts_blowup = 0;                    % re-initialize variables
            blockno = 1;
            lrate = lrate*DEFAULT_RESTART_FAC; % with lower learning rate
            weights = startweights;            % and original weight matrix
            oldweights = startweights;
            change = nochange;
            oldwtchange = zeros(chans,ncomps);
            delta=zeros(1,chans*ncomps);
            olddelta = delta;
            extblocks = urextblocks;
            prevweights = startweights;
            prevwtchange = zeros(chans,ncomps);
            lrates = zeros(1,maxsteps);
            bias = zeros(ncomps,1);

            signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1
            for k=1:nsub
                signs(k) = -1;
            end
            signs = diag(signs); % make a diagonal matrix
            oldsigns = zeros(size(signs));;

            if lrate> MIN_LRATE
                r = rank(data); % determine if data rank is too low 
                if r<ncomps
                    icaprintf(verb,fid,'Data has rank %d. Cannot compute %d components.\n',...
                        r,ncomps);
                    return
                else
                    icaprintf(verb,fid,...
                        'Lowering learning rate to %g and starting again.\n',lrate);
                end
            else
                icaprintf(verb,fid, ...
                    'runica(): QUITTING - weight matrix may not be invertible!\n');
                return;
            end
        else % if weights in bounds
            %
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step> 2
                angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
            end
            places = -floor(log10(nochange));
            icaprintf(verb,fid,'step %d - lrate %5f, wchange %8.8f, angledelta %4.1f deg\n', ...
                                step,      lrate,     change, degconst*angledelta);
            %
            %%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            changes = [changes change];
            oldweights = weights;
            %
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if degconst*angledelta > annealdeg,
                lrate = lrate*annealstep;          % anneal learning rate
                olddelta   = delta;                % accumulate angledelta until
                oldchange  = change;               %  annealdeg is reached
            elseif step == 1                     % on first step only
                olddelta   = delta;                % initialize
                oldchange  = change;
            end
            %
            %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if step >2 & change < nochange,      % apply stopping rule
                laststep=step;
                step=maxsteps;                  % stop when weights stabilize
            elseif change > DEFAULT_BLOWUP,      % if weights blow up,
                lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying
            end;                                 % with a smaller learning rate
        end; % end if weights in bounds

    end; % end while step < maxsteps (ICA Training) %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Compute ICA Weights
if biasflag & ~extended
    while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        timeperm=randperm(datalength); % shuffle data order at each step

        for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
            if strcmpi(interupt, 'on')
                drawnow;
                flag = getappdata(fig, 'run');
                if ~flag,
                    if ~isempty(fid), fclose(fid); end;
                    close; error('USER ABORT');
                end;
            end;
            
            u=weights*double(data(:,timeperm(t:t+block-1))) + bias*onesrow;
            y=1./(1+exp(-u));                                                
            weights = weights + lrate*(BI+(1-2*y)*u')*weights;               
            bias = bias + lrate*sum((1-2*y)')'; % for logistic nonlin. %
               
            if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                weights = weights + momentum*prevwtchange;
                prevwtchange = weights-prevweights;
                prevweights = weights;
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if max(max(abs(weights))) > MAX_WEIGHT
                wts_blowup = 1;
                change = nochange;
            end
            blockno = blockno + 1;
            if wts_blowup
                break
            end
        end % for t=1:block:lastt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~wts_blowup
            oldwtchange = weights-oldweights;
            step=step+1;
            %
            %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
            %
            lrates(1,step) = lrate;
            angledelta=0.;
            delta=reshape(oldwtchange,1,chans*ncomps);
            change=delta*delta';
        end
        %
        %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
        %
        if wts_blowup | isnan(change)|isinf(change),  % if weights blow up,
            icaprintf(verb,fid,'');
            step = 0;                          % start again
            change = nochange;
            wts_blowup = 0;                    % re-initialize variables
            blockno = 1;
            lrate = lrate*DEFAULT_RESTART_FAC; % with lower learning rate
            weights = startweights;            % and original weight matrix
            oldweights = startweights;
            change = nochange;
            oldwtchange = zeros(chans,ncomps);
            delta=zeros(1,chans*ncomps);
            olddelta = delta;
            extblocks = urextblocks;
            prevweights = startweights;
            prevwtchange = zeros(chans,ncomps);
            lrates = zeros(1,maxsteps);
            bias = zeros(ncomps,1);
            if lrate> MIN_LRATE
                r = rank(data); % determine if data rank is too low
                if r<ncomps
                    icaprintf(verb,fid,'Data has rank %d. Cannot compute %d components.\n',r,ncomps);
                    return
                else
                    icaprintf(verb,fid,'Lowering learning rate to %g and starting again.\n',lrate);
                end
            else
                icaprintf(verb,fid,'runica(): QUITTING - weight matrix may not be invertible!\n');
                return;
            end
        else % if weights in bounds
            %
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step> 2
                angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
            end
            places = -floor(log10(nochange));
            icaprintf(verb,fid,'step %d - lrate %5f, wchange %8.8f, angledelta %4.1f deg\n', ...
                                step,      lrate,     change, degconst*angledelta);
            %
            %%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            changes = [changes change];
            oldweights = weights;
            %
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if degconst*angledelta > annealdeg,
                lrate = lrate*annealstep;          % anneal learning rate
                olddelta   = delta;                % accumulate angledelta until
                oldchange  = change;               %  annealdeg is reached
            elseif step == 1                     % on first step only
                olddelta   = delta;                % initialize
                oldchange  = change;
            end
            %
            %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if step >2 & change < nochange,      % apply stopping rule
                laststep=step;
                step=maxsteps;                  % stop when weights stabilize
            elseif change > DEFAULT_BLOWUP,      % if weights blow up,
                lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying
            end;                                 % with a smaller learning rate
        end; % end if weights in bounds

    end; % end while step < maxsteps (ICA Training) %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Compute ICA Weights
if ~biasflag & extended
    while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        timeperm=randperm(datalength); % shuffle data order at each step through data

        for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
            if strcmpi(interupt, 'on')
                drawnow;
                flag = getappdata(fig, 'run');
                if ~flag,
                    if ~isempty(fid), fclose(fid); end;
                    close; error('USER ABORT');
                end;
            end;
            
            u=weights*double(data(:,timeperm(t:t+block-1))); % promote block to dbl
            y=tanh(u);                                                       %
            weights = weights + lrate*(BI-signs*y*u'-u*u')*weights;

            if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                weights = weights + momentum*prevwtchange;
                prevwtchange = weights-prevweights;
                prevweights = weights;
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if max(max(abs(weights))) > MAX_WEIGHT
                wts_blowup = 1;
                change = nochange;
            end
            if ~wts_blowup
                %
                %%%%%%%%%%% Extended-ICA kurtosis estimation %%%%%%%%%%%%%%%%%%%%%
                %while step < maxsteps
                if extblocks > 0 & rem(blockno,extblocks) == 0,
                    % recompute signs vector using kurtosis
                    if kurtsize < frames % 12-22-99 rand() size suggestion by M. Spratling
                        rp = fix(rand(1,kurtsize)*datalength);  % pick random subset
                        % Accout for the possibility of a 0 generation by rand
                        ou = find(rp == 0);
                        while ~isempty(ou) % 1-11-00 suggestion by J. Foucher
                            rp(ou) = fix(rand(1,length(ou))*datalength);
                            ou = find(rp == 0);
                        end
                        partact=weights*double(data(:,rp(1:kurtsize)));
                    else                                        % for small data sets,
                        partact=weights*double(data);           % use whole data
                    end
                    m2=mean(partact'.^2).^2;
                    m4= mean(partact'.^4);
                    kk= (m4./m2)-3.0;                           % kurtosis estimates
                    if extmomentum
                        kk = extmomentum*old_kk + (1.0-extmomentum)*kk; % use momentum
                        old_kk = kk;
                    end
                    signs=diag(sign(kk+signsbias));             % pick component signs
                    if signs == oldsigns,
                        signcount = signcount+1;
                    else
                        signcount = 0;
                    end
                    oldsigns = signs;
                    signcounts = [signcounts signcount];
                    if signcount >= SIGNCOUNT_THRESHOLD,
                        extblocks = fix(extblocks * SIGNCOUNT_STEP);% make kurt() estimation
                        signcount = 0;                             % less frequent if sign
                    end                                         % is not changing
                end % extblocks > 0 & . . .
            end % if ~wts_blowup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            blockno = blockno + 1;
            if wts_blowup
                break
            end
        end % for t=1:block:lastt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~wts_blowup
            oldwtchange = weights-oldweights;
            step=step+1;
            %
            %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
            %
            lrates(1,step) = lrate;
            angledelta=0.;
            delta=reshape(oldwtchange,1,chans*ncomps);
            change=delta*delta';
        end
        %
        %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
        %
        if wts_blowup | isnan(change)|isinf(change),  % if weights blow up,
            icaprintf(verb,fid,'');
            step = 0;                          % start again
            change = nochange;
            wts_blowup = 0;                    % re-initialize variables
            blockno = 1;
            lrate = lrate*DEFAULT_RESTART_FAC; % with lower learning rate
            weights = startweights;            % and original weight matrix
            oldweights = startweights;
            change = nochange;
            oldwtchange = zeros(chans,ncomps);
            delta=zeros(1,chans*ncomps);
            olddelta = delta;
            extblocks = urextblocks;
            prevweights = startweights;
            prevwtchange = zeros(chans,ncomps);
            lrates = zeros(1,maxsteps);
            bias = zeros(ncomps,1);
            signs = ones(1,ncomps);    % initialize signs to nsub -1, rest +1
            for k=1:nsub
                signs(k) = -1;
            end
            signs = diag(signs); % make a diagonal matrix
            oldsigns = zeros(size(signs));
            if lrate> MIN_LRATE
                r = rank(data); % find whether data rank is too low
                if r<ncomps
                    icaprintf(verb,fid,'Data has rank %d. Cannot compute %d components.\n',...
                        r,ncomps);
                    return
                else
                    icaprintf(verb,fid,...
                        'Lowering learning rate to %g and starting again.\n',lrate);
                end
            else
                icaprintf(verb,fid, ...
                    'runica(): QUITTING - weight matrix may not be invertible!\n');
                return;
            end
        else % if weights in bounds
            %
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step> 2
                angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
            end
            places = -floor(log10(nochange));
            icaprintf(verb,fid,'step %d - lrate %5f, wchange %8.8f, angledelta %4.1f deg\n', ...
                                step,      lrate,     change, degconst*angledelta);
            %
            %%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            changes = [changes change];
            oldweights = weights;
            %
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if degconst*angledelta > annealdeg,
                lrate = lrate*annealstep;          % anneal learning rate
                olddelta   = delta;                % accumulate angledelta until
                oldchange  = change;               %  annealdeg is reached
            elseif step == 1                     % on first step only
                olddelta   = delta;                % initialize
                oldchange  = change;
            end
            %
            %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if step >2 & change < nochange,      % apply stopping rule
                laststep=step;
                step=maxsteps;                  % stop when weights stabilize
            elseif change > DEFAULT_BLOWUP,      % if weights blow up,
                lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying
            end;                                 % with a smaller learning rate
        end; % end if weights in bounds

    end; % end while step < maxsteps (ICA Training) %%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute ICA Weights
if ~biasflag & ~extended
    while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        timeperm=randperm(datalength); % shuffle data order at each step

        for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%
            if strcmpi(interupt, 'on')
                drawnow;
                flag = getappdata(fig, 'run');
                if ~flag,
                    if ~isempty(fid), fclose(fid); end;
                    close; error('USER ABORT');
                end;
            end;
            u=weights*double(data(:,timeperm(t:t+block-1)));
            y=1./(1+exp(-u));                                                %
            weights = weights + lrate*(BI+(1-2*y)*u')*weights;

            if momentum > 0 %%%%%%%%% Add momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                weights = weights + momentum*prevwtchange;
                prevwtchange = weights-prevweights;
                prevweights = weights;
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if max(max(abs(weights))) > MAX_WEIGHT
                wts_blowup = 1;
                change = nochange;
            end
            
            blockno = blockno + 1;
            if wts_blowup
                break
            end
        end % for t=1:block:lastt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~wts_blowup
            oldwtchange = weights-oldweights;
            step=step+1;
            %
            %%%%%%% Compute and print weight and update angle changes %%%%%%%%%
            %
            lrates(1,step) = lrate;
            angledelta=0.;
            delta=reshape(oldwtchange,1,chans*ncomps);
            change=delta*delta';
        end
        %
        %%%%%%%%%%%%%%%%%%%%%% Restart if weights blow up %%%%%%%%%%%%%%%%%%%%
        %
        if wts_blowup | isnan(change)|isinf(change),  % if weights blow up,
            icaprintf(verb,fid,'');
            step = 0;                          % start again
            change = nochange;
            wts_blowup = 0;                    % re-initialize variables
            blockno = 1;
            lrate = lrate*DEFAULT_RESTART_FAC; % with lower learning rate
            weights = startweights;            % and original weight matrix
            oldweights = startweights;
            change = nochange;
            oldwtchange = zeros(chans,ncomps);
            delta=zeros(1,chans*ncomps);
            olddelta = delta;
            extblocks = urextblocks;
            prevweights = startweights;
            prevwtchange = zeros(chans,ncomps);
            lrates = zeros(1,maxsteps);
            bias = zeros(ncomps,1);
            
            if lrate> MIN_LRATE
                r = rank(data); % find whether data rank is too low
                if r<ncomps
                    icaprintf(verb,fid,'Data has rank %d. Cannot compute %d components.\n',...
                        r,ncomps);
                    return
                else
                    icaprintf(verb,fid,...
                        'Lowering learning rate to %g and starting again.\n',lrate);
                end
            else
                icaprintf(verb,fid, ...
                    'runica(): QUITTING - weight matrix may not be invertible!\n');
                return;
            end
        else % if weights in bounds
            %
            %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%%
            %
            if step> 2
                angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
            end
            places = -floor(log10(nochange));
            icaprintf(verb,fid,'step %d - lrate %5f, wchange %8.8f, angledelta %4.1f deg\n', ...
                                step,      lrate,     change, degconst*angledelta);
            %
            %%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            changes = [changes change];
            oldweights = weights;
            %
            %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if degconst*angledelta > annealdeg,
                lrate = lrate*annealstep;          % anneal learning rate
                olddelta   = delta;                % accumulate angledelta until
                oldchange  = change;               %  annealdeg is reached
            elseif step == 1                     % on first step only
                olddelta   = delta;                % initialize
                oldchange  = change;
            end
            %
            %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if step >2 & change < nochange,      % apply stopping rule
                laststep=step;
                step=maxsteps;                  % stop when weights stabilize
            elseif change > DEFAULT_BLOWUP,      % if weights blow up,
                lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying
            end;                                 % with a smaller learning rate
        end; % end if weights in bounds

    end; % end while step < maxsteps (ICA Training) %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Finalize Computed Data for Output
  
if strcmpi(interupt, 'on')
    close(fig);
end;


  if ~laststep
    laststep = step;
  end;
  lrates = lrates(1,1:laststep);           % truncate lrate history vector

  %
  %%%%%%%%%%%%%% Orient components towards max positive activation %%%%%%
  %
  if nargout > 6 | strcmp(posactflag,'on')
      % make activations from sphered and pca'd data; -sm 7/05
      % add back the row means removed from data before sphering
      if strcmp(pcaflag,'off')
          sr = sphere * rowmeans';
          for r = 1:ncomps
              data(r,:) = data(r,:)+sr(r); % add back row means 
          end
          data = weights*data; % OK in single
      else
          ser = sphere*eigenvectors(:,1:ncomps)'*rowmeans';
          for r = 1:ncomps
              data(r,:) = data(r,:)+ser(r); % add back row means 
          end
          data = weights*data; % OK in single
      end;
  end
  %
  % NOTE: Now 'data' are the component activations = weights*sphere*raw_data
  %
  %
  %%%%%%%%%%%%%% If pcaflag, compose PCA and ICA matrices %%%%%%%%%%%%%%%
  %
  if strcmp(pcaflag,'on')
    icaprintf(verb,fid,'Composing the eigenvector, weights, and sphere matrices\n');
    icaprintf(verb,fid,'  into a single rectangular weights matrix; sphere=eye(%d)\n'...
                                                                  ,chans);
    weights= weights*sphere*eigenvectors(:,1:ncomps)'; 
    sphere = eye(urchans);
  end
  %
  %%%%%% Sort components in descending order of max projected variance %%%%
  %
  icaprintf(verb,fid,'Sorting components in descending order of mean projected variance ...\n');
  %
  %%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % meanvar  = zeros(ncomps,1);      % size of the projections
  if ncomps == urchans % if weights are square . . .
      winv = inv(weights*sphere);
  else
      icaprintf(verb,fid,'Using pseudo-inverse of weight matrix to rank order component projections.\n');
      winv = pinv(weights*sphere);
  end
  %
  % compute variances without backprojecting to save time and memory -sm 7/05
  %
  meanvar = sum(winv.^2).*sum((data').^2)/((chans*frames)-1); % from Rey Ramirez 8/07
  %
  %%%%%%%%%%%%%% Sort components by mean variance %%%%%%%%%%%%%%%%%%%%%%%%
  %
  [sortvar, windex] = sort(meanvar);
  windex = windex(ncomps:-1:1); % order large to small 
  meanvar = meanvar(windex);
  %
  %%%%%%%%%%%% re-orient max(abs(activations)) to >=0 ('posact') %%%%%%%%
  %
  if strcmp(posactflag,'on') % default is now off to save processing and memory
      icaprintf(verb,fid,'Making the max(abs(activations)) positive ...\n');
      [tmp ix] = max(abs(data')); % = max abs activations
      signsflipped = 0;
      for r=1:ncomps
         if sign(data(r,ix(r))) < 0
            if nargout>6  % if activations are to be returned (only)
               data(r,:) = -1*data(r,:);  % flip activations so max(abs()) is >= 0
            end
            winv(:,r) = -1*winv(:,r);  % flip component maps
            signsflipped = 1;
         end
      end
      if signsflipped == 1
          weights = pinv(winv)*inv(sphere); % re-invert the component maps
      end
       
      % [data,winvout,weights] = posact(data,weights); % overwrite data with activations
      % changes signs of activations (now = data) and weights 
      % to make activations (data) net rms-positive
      % can call this outside of runica() - though it is inefficient!
  end
  % 
  %%%%%%%%%%%%%%%%%%%%% Filter data using final weights %%%%%%%%%%%%%%%%%%
  %
  if nargout>6, % if activations are to be returned
      icaprintf(verb,fid,'Permuting the activation wave forms ...\n');
      data = data(windex,:); % data is now activations -sm 7/05
  else
      clear data
  end
  weights = weights(windex,:);% reorder the weight matrix
  bias  = bias(windex);       % reorder them
  signs = diag(signs);        % vectorize the signs matrix
  signs = signs(windex);      % reorder them
  
  if ~isempty(fid), fclose(fid); end; % close logfile

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    
return

% printing functions
% ------------------
function icaprintf(verb,fid, varargin);
    if verb
        if ~isempty(fid)
            fprintf(fid, varargin{:});
        end;        
        fprintf(varargin{:});
    end;
