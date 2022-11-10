function [data] = ft_timelocksimulation(cfg)

% FT_TIMELOCKSIMULATION computes simulated data that consists of multiple trials in
% with each trial contains an event-related potential or field. Following
% construction of the time-locked signal in each trial by this function, the signals
% can be passed into FT_TIMELOCKANALYSIS to obtain the average and the variance.
%
% Use as
%   [data] = ft_timelockstatistics(cfg)
% which will return a raw data structure that resembles the output of
% FT_PREPROCESSING.
%
% The number of trials and the time axes of the trials can be specified by
%   cfg.fsample    = simulated sample frequency (default = 1000)
%   cfg.trllen     = length of simulated trials in seconds (default = 1)
%   cfg.numtrl     = number of simulated trials (default = 10)
%   cfg.baseline   = number (default = 0.3)
% or by
%   cfg.time       = cell-array with one time axis per trial, which are for example obtained from an existing dataset
%
% The signal is constructed from three underlying functions. The shape is
% controlled with
%   cfg.s1.numcycli = number (default = 1)
%   cfg.s1.ampl     = number (default = 1.0)
%   cfg.s2.numcycli = number (default = 2)
%   cfg.s2.ampl     = number (default = 0.7)
%   cfg.s3.numcycli = number (default = 4)
%   cfg.s3.ampl     = number (default = 0.2)
%   cfg.noise.ampl  = number (default = 0.1)
% Specifying numcycli=1 results in a monophasic signal, numcycli=2 is a biphasic,
% etc. The three signals are scaled to the indicated amplitude, summed up and a
% certain amount of noise is added.
%
% Other configuration options include
%   cfg.numchan     = number (default = 5)
%
% See also FT_TIMELOCKANALYSIS, FT_TIMELOCKSTATISTICS, FT_FREQSIMULATION,
% FT_DIPOLESIMULATION, FT_CONNECTIVITYSIMULATION

% Copyright (C) 2016-2020, Robert Oostenveld
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

% Ideas to extend this function
% - add some jitter to each signal on each trial
% - add some amplitude variation on each signal on each trial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% get the options
cfg.numchan   = ft_getopt(cfg, 'numchan', 5);
cfg.time      = ft_getopt(cfg, 'time', []);
if isempty(cfg.time)
  cfg.fsample   = ft_getopt(cfg, 'fsample', 1000);
  cfg.trllen    = ft_getopt(cfg, 'trllen', 1);
  cfg.numtrl    = ft_getopt(cfg, 'numtrl', 10);
  cfg.baseline  = ft_getopt(cfg, 'baseline', 0.3);
else
  cfg.trllen    = length(cfg.time{1});        % must be identical for all trials
  cfg.fsample   = 1/mean(diff(cfg.time{1}));  % determine from time-axis
  cfg.numtrl    = length(cfg.time);
end

cfg.s1 = ft_getopt(cfg, 's1');
cfg.s2 = ft_getopt(cfg, 's2');
cfg.s3 = ft_getopt(cfg, 's3');
cfg.noise = ft_getopt(cfg, 'noise');

cfg.s1.numcycli = ft_getopt(cfg.s1, 'numcycli', 1);
cfg.s1.ampl     = ft_getopt(cfg.s1, 'ampl', 1.0);
cfg.s2.numcycli = ft_getopt(cfg.s2, 'numcycli', 2);
cfg.s2.ampl     = ft_getopt(cfg.s2, 'ampl', 0.7);
cfg.s3.numcycli = ft_getopt(cfg.s3, 'numcycli', 4);
cfg.s3.ampl     = ft_getopt(cfg.s3, 'ampl', 0.2);
cfg.noise.ampl  = ft_getopt(cfg.noise, 'ampl', 0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the simulated timeseries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(cfg.time)
  % use the user-supplied time vectors
  timevec = cfg.time;
else
  % give the user some feedback
  ft_debug('using %f as samping frequency', cfg.fsample);
  ft_debug('using %d trials of %f seconds long', cfg.numtrl, cfg.trllen);
  nsample = round(cfg.trllen*cfg.fsample);
  timevec = cell(1, cfg.numtrl);
  for iTr = 1:cfg.numtrl
    timevec{iTr} = (((1:nsample)-1)/cfg.fsample) - cfg.baseline;
  end
end

data = [];
data.label = {};
for i=1:cfg.numchan
  data.label{i} = sprintf('%d', i);
end
data.time  = timevec;
data.trial = cell(1, numel(data.time));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct each of the trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:numel(data.time)
  nsample   = length(data.time{i});
  nbaseline = sum(data.time{i}<0);
  sel       = (nbaseline+1):nsample;
  
  % start with a prototype monophasic signal
  % it will become multiphasic by taking the n-th derivative
  signal = gausswin(nsample-nbaseline, 4)';
  signal = signal .* hanning(nsample-nbaseline)';
  signal = signal-signal(1);
  
  if ~isempty(cfg.s1)
    signal1 = signal;
    countdown = cfg.s1.numcycli;
    while countdown>1
      signal1 = gradient(signal1);
      signal1 = signal1-signal1(1);
      countdown = countdown - 1;
    end
    signal1 = signal1/max(abs(signal1)) * cfg.s1.ampl;
    % figure; plot(signal1); title('signal 1')
  end
  
  if ~isempty(cfg.s2)
    signal2 = signal;
    countdown = cfg.s2.numcycli;
    while countdown>1
      signal2 = gradient(signal2);
      signal2 = signal2-signal2(1);
      countdown = countdown - 1;
    end
    signal2 = signal2/max(abs(signal2)) * cfg.s2.ampl;
    % figure; plot(signal2); title('signal 2')
  end
  
  if ~isempty(cfg.s3)
    signal3 = signal;
    countdown = cfg.s3.numcycli;
    while countdown>1
      signal3 = gradient(signal3);
      signal3 = signal3-signal3(1);
      countdown = countdown - 1;
    end
    signal3 = signal3/max(abs(signal3)) * cfg.s3.ampl;
    % figure; plot(signal3); title('signal 3')
  end
  
  % start with an empty data matrix for this trial
  dat = zeros(numel(data.label), nsample);
  
  for j=1:cfg.numchan
    if ~isempty(cfg.s1)
      dat(j,sel) = dat(j,sel) + signal1;
    end
    if ~isempty(cfg.s2)
      dat(j,sel) = dat(j,sel) + signal2;
    end
    if ~isempty(cfg.s3)
      dat(j,sel) = dat(j,sel) + signal3;
    end
    if ~isempty(cfg.noise)
      % the signal is the same, but the noise is different for each channel
      dat(j,:) = dat(j,:) + randn(1,nsample)*cfg.noise.ampl;
    end
  end % for numchan
  
  data.trial{i} = dat;
  data.time{i} = timevec{i};
end % for each trial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the bookkeeping at the end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_postamble debug
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
