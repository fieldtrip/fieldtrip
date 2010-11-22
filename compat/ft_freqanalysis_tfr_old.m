function [freq] = ft_freqanalysis_tfr(cfg, data);

% FT_FREQANALYSIS_TFR computes time-frequency representations of single-trial
% data using a convolution in the time-domain with Morlet's wavelets.
%
% Use as
%   [freq] = ft_freqanalysis(cfg, data)
%
% The data should be organised in a structure as obtained from
% the FT_PREPROCESSING function. The configuration should be according to
%   cfg.method        = method used for frequency or time-frequency decomposition
%                       see FT_FREQANALYSIS for details
%   cfg.foi           = vector 1 x numfoi, frequencies of interest
%   cfg.waveletwidth  = 'width' of wavelets expressed in cycles (default = 7)
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                       see FT_CHANNELSELECTION for details
%   cfg.downsample    = ratio for downsampling, which occurs after convolution (default = 1)
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.keeptrials    = 'yes' or 'no', return individual trials or average (default = 'no')
%
% See also FT_FREQANALYSIS

% Undocumented local options:
% cfg.latency
% cfg.output

% Copyright (C) 2003, Ole Jensen, FCDC
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
% $Id: ft_freqanalysis_tfr.m 1980 2010-10-27 10:45:10Z jansch $

fieldtripdefs

% ensure that this function is started as a subfunction of the FT_FREQANALYSIS wrapper
if ~exist('OCTAVE_VERSION')
  [s, i] = dbstack;
  if length(s)>1
    [caller_path, caller_name, caller_ext] = fileparts(s(2).name);
  else
    caller_path = '';
    caller_name = '';
    caller_ext  = '';
  end
  % evalin('caller', 'mfilename') does not work for Matlab 6.1 and 6.5
  if ~strcmp(caller_name, 'ft_freqanalysis')
    error(['you should call FREQANALYSIS, instead of the ' upper(mfilename) ' subfunction']);
  end
end

% set the defaults
if ~isfield(cfg, 'method'),         cfg.method  = 'tfr';          end
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';          end
if ~isfield(cfg, 'latency'),        cfg.latency = 'minperlength'; end
if ~isfield(cfg, 'keeptrials'),     cfg.keeptrials   = 'no';      end
if ~isfield(cfg, 'waveletwidth'),   cfg.waveletwidth = 7;         end
if ~isfield(cfg, 'downsample'),     cfg.downsample   = 1;         end
if ~isfield(cfg, 'feedback'),       cfg.feedback     = 'text';    end

if isfield(cfg, 'output') && strcmp(cfg.output, 'powandcsd'),
  error('This function does not compute cross-spectra\n');
end

% determine the channels of interest
cfg.channel = ft_channelselection(cfg.channel, data.label);
chansel     = match_str(data.label, cfg.channel);

% determine the duration of each trial
ntrial = length(data.trial);
nchan = size(data.trial{1}, 1);

for i=1:ntrial
  nsampl(i)         = size(data.trial{i}, 2);
  begsamplatency(i) = min(data.time{i});
  endsamplatency(i) = max(data.time{i});
end;

if cfg.downsample > 1
  % perform a decimation of the input data
  warning('decimating the input data, better is to use RESAMPLEDATA');
  for k=1:ntrial
    dTmp = data.trial{k};
    data.trial{k} = dTmp(:,1:cfg.downsample:end);
    tTmp = data.time{k};
    data.time{k} = tTmp(1:cfg.downsample:end);
  end
  data.fsample = data.fsample / cfg.downsample;
end

% automatically determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];

% latency window for averaging and variance computation is given in seconds
if (strcmp(cfg.latency, 'minperlength'))
  cfg.latency = [];
  cfg.latency(1) = minperlength(1);
  cfg.latency(2) = minperlength(2);
end
if (strcmp(cfg.latency, 'prestim'))
  cfg.latency = [];
  cfg.latency(1) = minperlength(1);
  cfg.latency(2) = 0;
end
if (strcmp(cfg.latency, 'poststim'))
  cfg.latency = [];
  cfg.latency(1) = 0;
  cfg.latency(2) = minperlength(2);
end

M = waveletfam(cfg.foi,data.fsample,cfg.waveletwidth);

ft_progress('init', cfg.feedback, 'convolving wavelets');

for i=1:ntrial
  indicvect = data.time{i};
  ft_progress(i/ntrial, 'convolving wavelets, trial %d of %d\n', i, ntrial);

  %for average and variance
  begsampl = nearest(indicvect,cfg.latency(1));
  endsampl = nearest(indicvect,cfg.latency(2));

  numsamples(i) = endsampl-begsampl+1;

  if (i==1)
    % allocate memory to hold the resulting powerspectra
    if strcmp(cfg.keeptrials, 'yes')
      freq.powspctrm = zeros(ntrial,nchan,length(cfg.foi),ceil((endsampl-begsampl+1)/cfg.downsample));
    else
      freq.powspctrm = zeros(nchan,length(cfg.foi),ceil((endsampl-begsampl+1)/cfg.downsample));
    end
  end;

  dat = data.trial{i}(chansel,begsampl:endsampl);
  for k=1:size(dat,1)
    for j=1:length(cfg.foi)
      cTmp = conv(dat(k,:),M{j});
      cTmp = 2*(abs(cTmp).^2)/data.fsample;
      cTmp = cTmp(ceil(length(M{j})/2):length(cTmp)-floor(length(M{j})/2));
      cTmp = cTmp(:,1:cfg.downsample:end);
      if strcmp(cfg.keeptrials, 'yes')
        freq.powspctrm(i,k,j,:) = cTmp';
      else
        freq.powspctrm(k,j,:) = squeeze(freq.powspctrm(k,j,:)) + cTmp';  % compute the running sum
      end
    end
  end

end %for ntrial

ft_progress('close');

if strcmp(cfg.keeptrials, 'yes')
  freq.dimord    = 'rpt_chan_freq_time';
else
  freq.dimord    = 'chan_freq_time';
  freq.powspctrm = freq.powspctrm / ntrial;  % compute the average
end
freq.label     = cfg.channel;
freq.freq      = cfg.foi;
freq.time      = indicvect(1:cfg.downsample:end);

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: ft_freqanalysis_tfr.m 1980 2010-10-27 10:45:10Z jansch $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
freq.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for waveletanalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = waveletfam(foi,Fs,width)
dt = 1/Fs;
for k=1:length(foi)
  sf = foi(k)/width;
  st = 1/(2*pi*sf);
  toi=-3.5*st:dt:3.5*st;
  A = 1/sqrt(st*sqrt(pi));
  M{k}= A*exp(-toi.^2/(2*st^2)).*exp(i*2*pi*foi(k).*toi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for waveletanalysis
%
% Return a vector containing the energy as a
% function of time for frequency f. The energy
% is calculated using Morlet's wavelets.
% s : signal
% Fs: sampling frequency
% width: width of Morlet wavelet (>= 5 suggested).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = energyvec(f,s,Fs,width)
dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);
t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);
size(m)
y = conv(s,m);
y = (2*abs(y)/Fs).^2;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction for waveletanalysis
%
% Morlet's wavelet for frequency f and time t.
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet.
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY
%
% Ole Jensen, August 1998
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = morlet(f,t,width)
sf = f/width;
st = 1/(2*pi*sf);
A = 1/sqrt(st*sqrt(pi));
y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);
