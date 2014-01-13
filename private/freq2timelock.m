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
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
if isfield(freq, 'cfg'), timelock.cfg = freq.cfg; end
timelock.dimord = 'chan_time';

if isfield(freq, 'grad')
  timelock.grad = freq.grad;
end

if isfield(freq, 'elec')
  timelock.elec = freq.elec;
end

