function [data] = ft_connectivitysimulation(cfg)

% FT_CONNECTIVITYSIMULATION
%
% Use as
%   data = ft_connectivitysimulation(cfg)
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
% See also FT_FREQSIMULATION, FT_DIPOLESIMULATION, FT_SPIKESIMULATION

% Copyright (C) 2009, Donders Institute for Brain, Cognition and Behaviour

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

  fltpad = 50; %hard coded to avoid filtering artifacts
  delay  = cfg.delay;
  delay  = delay - min(delay(:)); %make explicitly >= 0
  maxdelay = max(delay(:));

  if iscell(cfg.mix),
    %each trial has different mix
    mix = cfg.mix;            
  else
    %make cell-array out of mix
    tmpmix = cfg.mix;
    mix    = cell(1,cfg.ntrials);
    for tr = 1:cfg.ntrials
      mix{1,tr} = tmpmix;
    end
  end
  
  nmixsignal = size(mix{1}, 2); %number of "mixing signals"
  nsignal    = size(mix{1}, 1);

  if numel(size(mix{1}))==2,
    %mix is static, no function of time
    for tr = 1:cfg.ntrials
      mix{tr} = mix{tr}(:,:,ones(1,nsmp+maxdelay));
    end
  elseif numel(size(mix{1}))==3 && size(mix{1},3)==nsmp,
    %mix changes with time
    for tr = 1:cfg.ntrials
      mix{tr} = cat(3,mix{tr},mix{tr}(:,:,nsmp*ones(1,maxdelay))); 
    end
    %FIXME think about this
    %due to the delay the mix cannot be defined instantaneously with respect to all signals
  end
    
  for tr = 1:cfg.ntrials
    mixsignal = randn(nmixsignal,  nsmp + 2*fltpad + maxdelay);
    mixsignal = preproc(mixsignal, label, cfg.fsample, cfg, -fltpad, fltpad, fltpad);
    tmp       = zeros(cfg.nsignal, nsmp);
    for i=1:cfg.nsignal
      for j=1:nmixsignal
        begsmp   = 1    + delay(i,j);
        endsmp   = nsmp + delay(i,j);
        tmpmix   = reshape(mix{tr}(i,j,:),[1 nsmp+maxdelay]) .* mixsignal(j,:);
        tmp(i,:) = tmp(i,:) + tmpmix(begsmp:endsmp);
      end
    end
    trial{tr} = tmp;
    
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
cfg.version.id   = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
data.cfg     = cfg;
