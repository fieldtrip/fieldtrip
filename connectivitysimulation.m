function [data] = connectivitysimulation(cfg)

% CONNECTIVITYSIMULATION
%
% Use as
%   data = connectivitysimulation(cfg)
%
% cfg.method = string, can be xxx
%
% cfg.ntrials     = scalar
% cfg.triallength = in seconds
% cfg.fsample     = in Hz
% 
% cfg.bpfilter = 'yes' (or 'no')
% cfg.bpfreq   = [bplow bphigh]
% 
% cfg.nsignal     = scalar, number of signals
% cfg.covmat      = covariance matrix between the signals
% cfg.delay       = delay vector between the signals in samples
%
% See also FREQSIMULATION, DIPOLESIMULATION, SPIKESIMULATION

% Copyright (C) 2009, Donders Institute for Brain, Cognition and Behaviour
%
% $Log: connectivitysimulation.m,v $
% Revision 1.11  2009/10/15 11:33:47  jansch
% added noisecov for ar method
%
% Revision 1.10  2009/10/12 13:23:10  andbas
% took out the bandpass etc for ar method
%
% Revision 1.9  2009/10/09 15:17:32  andbas
% Added bandpass filtering and additional error terms to the AR method estimate.
% Does not appear to be working...
%
% Revision 1.8  2009/10/08 07:03:52  jansch
% changed initialization for ar-method. added scaling factor for noise
%
% Revision 1.7  2009/10/08 06:59:52  jansch
% removed some duplicate code
%
% Revision 1.6  2009/10/07 18:58:17  jansch
% added not yet functional ar method. updated linear_mix method, implemented
% by andre
%
% Revision 1.5  2009/10/06 16:21:17  andbas
% *** empty log message ***
%
% Revision 1.4  2009/10/02 13:50:31  andbas
% replaced calls to blc with preproc_baselinecorrect
%
% Revision 1.3  2009/09/30 12:49:03  jansch
% first working implementation
%
% Revision 1.2  2009/09/28 15:39:47  jansch
% first working implementation (simple)
%
% Revision 1.1  2009/09/28 11:22:50  roboos
% created initial framework, sofar only some documentation
% the idea is that Andre and Jan-Mathijs start filling it with code
%

% check input configuration for the generally applicable options
cfg = checkconfig(cfg, 'required', {'nsignal' 'ntrials' 'triallength' 'fsample' 'method'});

% method specific defaults
switch cfg.method
case {'linear_mix'}
  %method specific defaults
  if ~isfield(cfg, 'bpfilter'), cfg.bpfilter = 'yes';   end
  if ~isfield(cfg, 'bpfreq'),   cfg.bpfreq   = [15 25]; end
  if ~isfield(cfg, 'blc'),      cfg.blc      = 'yes';   end
  if ~isfield(cfg, 'absnoise'), cfg.absnoise = 1;       end
  cfg = checkconfig(cfg, 'required', {'mix' 'delay'});
case {'mvnrnd'}
  if ~isfield(cfg, 'bpfilter'), cfg.bpfilter = 'yes';   end
  if ~isfield(cfg, 'bpfreq'),   cfg.bpfreq   = [15 25]; end
  if ~isfield(cfg, 'blc'),      cfg.blc      = 'yes';   end
  if ~isfield(cfg, 'absnoise'), cfg.absnoise = 1;       end
  cfg = checkconfig(cfg, 'required', {'covmat' 'delay'}); 
case {'ar'}
  cfg = checkconfig(cfg, 'required', {'params' 'noisecov'});
otherwise
end

trial = cell(1, cfg.ntrials);
time  = cell(1, cfg.ntrials);
nsmp  = round(cfg.triallength*cfg.fsample);
tim   = [0:nsmp-1]./cfg.fsample;

% create the labels
for k = 1:cfg.nsignal
  label{k,1} = ['signal',num2str(k,'%03d')];
end

switch cfg.method
case {'linear_mix'}
  fltpad = 100; %hard coded
  delay = cfg.delay;
  mix = cfg.mix;            
  nmixsignal  = size(cfg.mix,2); %number of "mixing signals"
  delay = delay - min(delay(:));
  for tr = 1:cfg.ntrials
    mixsignal = randn(nmixsignal, nsmp + 2*fltpad + max(delay(:)));
    tmp   = zeros(cfg.nsignal, nsmp+2*fltpad);
    for i=1:cfg.nsignal
      for j=1:nmixsignal
      begsmp = 1        + delay(i,j);
      endsmp = nsmp     + delay(i,j) + 2*fltpad;
      tmp(i,:) = tmp(i,:) + mix(i,j) * mixsignal(j,begsmp:endsmp);
      end
    end
    % apply preproc
    newtmp = preproc(tmp, label, cfg.fsample, cfg, -fltpad, fltpad, fltpad);
    trial{tr} = newtmp;
    % add some noise
    trial{tr} = preproc_baselinecorrect(trial{tr} + cfg.absnoise*randn(size(trial{tr})));
    % define time axis for this trial
    time{tr}  = tim;
  end
  
case {'mvnrnd'}
  fltpad = 100; %hard coded
  
  shift = max(cfg.delay(:,1)) - cfg.delay(:,1);
  for k = 1:cfg.ntrials
    % create the multivariate time series plus some padding
    tmp = mvnrnd(zeros(1,cfg.nsignal), cfg.covmat, nsmp+2*fltpad+max(shift))';
  
    % add the delays
    newtmp = zeros(cfg.nsignal, nsmp+2*fltpad);
    for kk = 1:cfg.nsignal
      begsmp =      + shift(kk) + 1;
      endsmp = nsmp + 2*fltpad + shift(kk);
      newtmp(kk,:) = preproc_baselinecorrect(tmp(kk,begsmp:endsmp));
    end
  
    % apply preproc
    newtmp = preproc(newtmp, label, cfg.fsample, cfg, -fltpad, fltpad, fltpad);
  
    trial{k} = newtmp;
    
    % add some noise
    trial{k} = preproc_baselinecorrect(trial{k} + cfg.absnoise*randn(size(trial{k})));
  
    % define time axis for this trial
    time{k}  = tim;
  end
case {'ar'}
  nlag    = size(cfg.params,3);
  nsignal = cfg.nsignal;
  params  = zeros(nlag*nsignal, nsignal);
  for k = 1:nlag
    params(((k-1)*nsignal+1):k*nsignal,:) = cfg.params(:,:,k);
  end 
  for k = 1:cfg.ntrials
    tmp   = zeros(nsignal, nsmp+nlag);
    noise  = mvnrnd(zeros(nsignal,1), cfg.noisecov, nsmp+nlag)';
    state0 = zeros(nsignal*nlag, 1);
    for m = 1:nlag
      indx = ((m-1)*nsignal+1):m*nsignal;
      state0(indx) = params(indx,:)'*noise(:,m);    
    end
    tmp(:,1:nlag) = fliplr(reshape(state0, [nsignal nlag]));  
    
    for m = (nlag+1):(nsmp+nlag)
       state0    = reshape(fliplr(tmp(:,(m-nlag):(m-1))), [nlag*nsignal 1]);
       tmp(:, m) = params'*state0 + noise(:,m); 
    end
    trial{k} = tmp(:,nlag+1:end);
    time{k}  = tim;
  end  

otherwise
  error('unknown method');
end

% create the output data
data         = [];
data.trial   = trial;
data.time    = time;
data.fsample = cfg.fsample;
data.label   = label;

% add version details to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: connectivitysimulation.m,v 1.11 2009/10/15 11:33:47 jansch Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
data.cfg     = cfg;
