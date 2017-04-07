function data = ft_spikesimulation(cfg)

% FT_SPIKESIMULATION simulates a spiketrain with a structures timing of the
% neuronal firing.
%
% Use as
%   data = ft_spikesimulation(cfg)
% and please look in the code for the cfg details.
%
% See also FT_DIPOLESIMULATION, FT_FREQSIMULATION

% Copyright (C) 2007, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig

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
    lfp(i,:) = ft_preproc_bandpassfilter(randn(1,nsample), fsample, cfg.bpfreq);
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

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble provenance data
ft_postamble history    data

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

