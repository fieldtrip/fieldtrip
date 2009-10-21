function [timelock, cfg] = freq2timelock(cfg, freq);

% FREQ2TIMELOCK  transform the frequency data into something
% on which the timelocked source reconstruction methods can
% perform their trick.
%
% This function also performs frequency and channel selection, using
%   cfg.frequency
%   cfg.channel
%
% After source reconstruction, you should use TIMELOCK2FREQ.

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: freq2timelock.m,v $
% Revision 1.5  2009/02/02 13:22:27  jansch
% changed 'chancmb' into 'chan' (line 32)
%
% Revision 1.4  2006/05/10 08:19:45  roboos
% added dimord to the output
%
% Revision 1.3  2006/02/23 10:28:17  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.2  2006/02/01 12:26:04  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.1  2005/10/14 15:50:08  roboos
% new implementation, used by dipolefitting in case of frequency or ICA data
%

if isfield(freq, 'fourierspctrm')
  fprintf('constructing real/imag data representation from single trial fourier representation\n');
  % select the complex amplitude at the frequency of interest
  cdim = dimnum(freq.dimord, 'chan');  % should be 2
  fdim = dimnum(freq.dimord, 'freq');     % should be 3
  fbin = nearest(freq.freq, cfg.frequency);
  cfg.frequency = freq.freq(fbin);
  if cdim==2 && fdim==3
    % other dimords are not supported, since they do not occur
    spctrm = dimindex(freq.fourierspctrm, fdim, fbin)';
  end
  % select the desired channels in the data
  cfg.channel = channelselection(cfg.channel, freq.label);
  [dum, chansel] = match_str(cfg.channel, freq.label);
  spctrm = spctrm(chansel,:);
  % concatenate the real and imaginary part
  avg = [real(spctrm) imag(spctrm)];
elseif isfield(freq, 'crsspctrm')
  fprintf('constructing real/imag data representation from csd matrix\n');
  % hmmm... I have no idea whether this is correct
  cfg.channel = channelselection(cfg.channel, freq.label);
  % this subfunction also takes care of the channel selection
  [Cf, Cr, Pr, Ntrials, dum] = prepare_freq_matrices(cfg, freq);
  cfg.frequency = dum.frequency;
  cfg.channel   = dum.channel;
  if length(size(Cf))==3
    % average the cross-spectrum over trials
    Cf = squeeze(mean(Cf,1));
  end
  % reconstruct something that spans the same space as the fft of the data, hmmm...
  [u, s, v] = svd(Cf);
  spctrm = u * sqrt(s);
  % concatenate the real and imaginary part
  avg = [real(spctrm) imag(spctrm)];
else
  error('unknown representation of frequency domain data');
end

timelock        = [];
timelock.avg    = avg;
timelock.label  = cfg.channel;
timelock.time   = 1:size(timelock.avg,2);
timelock.cfg    = freq.cfg;
timelock.dimord = 'chan_time';

if isfield(freq, 'grad')
  timelock.grad = freq.grad;
end

if isfield(freq, 'elec')
  timelock.elec = freq.elec;
end

