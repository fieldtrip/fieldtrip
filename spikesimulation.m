function data = spikesimulation(cfg)

% SPIKESIMULATION
%
% Use as
%   data = spikesimulation(cfg)
% and please look in the code for the cfg details.
%
% See also DIPOLESIMULATION, FREQSIMULATION

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: spikesimulation.m,v $
% Revision 1.7  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.6  2008/06/17 16:13:18  sashae
% now using preproc_modules
%
% Revision 1.5  2007/12/12 10:03:31  roboos
% fixed a cvs conflict, not sure whether this version is operational
%
% Revision 1.4  2007/11/05 09:43:54  roboos
% some cleanup
%
% Revision 1.3  2007/11/01 16:33:57  roboos
% added cfg to output
%
% Revision 1.2  2007/11/01 16:32:21  roboos
% updated documentation, some small fixes
%

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'trlduration'),  cfg.trlduration = 1; end %  in seconds
if ~isfield(cfg, 'nlfpchan'),     cfg.nlfpchan    = 10; end
if ~isfield(cfg, 'nspikechan'),   cfg.nspikechan  = 10; end
if ~isfield(cfg, 'spikerate'),    cfg.spikerate   = ones(cfg.nspikechan,1)*10; end
if ~isfield(cfg, 'ntrial'),       cfg.ntrial      = 10; end
if ~isfield(cfg, 'bpfreq'),       cfg.bpfreq      = [40 80]; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cfg.spikerate   = [ 10 50 100 200];
% cfg.nspikechan  = 4;
% nlfpchan        = 2;
% cfg.ntrial      = 10;
% cfg.bpfreq      = [40 80];
% cfg.trlduration = 1;

spikemix        = eye(cfg.nspikechan, cfg.nlfpchan);
spikemix(5,1:2) = [0.5 0.5]; % mix with lfp chan 1 and 2
spikemix(6,1:2) = [0.5 0.5]; % mix with lfp chan 1 and 2

% cfg.spikerate  = 100*ones(1,cfg.nspikechan);  % for each channel, per second
% cfg.spikerate  = [100 100 300 300 100 300];
fsample    = 1000;
nsample    = cfg.trlduration*fsample;
nchan      = cfg.nlfpchan + cfg.nspikechan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start with empty data
data = [];
data.fsample = fsample;

% create channel labels
k = 1;
for i=1:cfg.nlfpchan,    data.label{k} = sprintf('lfp%02d',   i); k=k+1; end
for i=1:cfg.nspikechan,  data.label{k} = sprintf('spike%02d', i); k=k+1; end
data.label = data.label(:);

for t=1:cfg.ntrial
  fprintf('creating trial %d\n', t);
  data.time{t} = (1:nsample)./fsample;
  lfp   = zeros(cfg.nlfpchan, nsample);
  spike = zeros(cfg.nspikechan, nsample);
  
  for i=1:cfg.nlfpchan, 
    lfp(i,:) = preproc_bandpassfilter(randn(1,nsample), fsample, cfg.bpfreq);
  end

  for i=1:cfg.nspikechan
    % the spikes are generated from a probabilistic mix of the LFP channels
    x = spikemix(i,:) * lfp;
    % apply a non-linear mapping to the spike probablility, output should be defined between 0 and 1
    x = mapping(x);
    % normalize the total probability to one
    x = x./sum(x);
    % randomly assign the spikes over the time series
    spike(i,:) = ((cfg.spikerate(i)*nsample/fsample)*x)>=rand(size(x));
  end

  data.time{t}  = (1:nsample)./fsample;
  data.trial{t} = [lfp; spike];
  clear lfp spike
end

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: spikesimulation.m,v 1.7 2008/09/22 20:17:44 roboos Exp $';
% remember the configuration details
data.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to strengthen the temporal strucure in the spike train
% by a non-linear modulation of the spike probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = mapping(x)
x = 6*(x-mean(x))/std(x);
y = 1./(1+exp(-x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IIRFILTER simple brick wall filter in the frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = myiirfilter(s, fs, fb, a)
f  = fft(s);
ax = linspace(0, fs, length(s));
fl = nearest(ax, fb(1))-1;
fh = nearest(ax, fb(2))+1;
f(1:fl)   = a.*f(1:fl);
f(fh:end) = a.*f(fh:end);
s = real(ifft(f));

