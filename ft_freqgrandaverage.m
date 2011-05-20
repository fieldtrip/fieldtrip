function [grandavg] = ft_freqgrandaverage(cfg, varargin);

% FT_FREQGRANDAVERAGE computes the average powerspectrum or time-frequency spectrum
% over multiple subjects
%
% Use as
%   [grandavg] = ft_freqgrandaverage(cfg, freq1, freq2, freq3...)
%
% The input data freq1..N are obtained from either FT_FREQANALYSIS with
% keeptrials=no or from FT_FREQDESCRIPTIVES. The configuration structure
% can contain
%  cfg.keepindividual = 'yes' or 'no' (default = 'no')
%  cfg.foilim         = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%  cfg.toilim         = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
%  cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                       see FT_CHANNELSELECTION for details
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure. For this particular function, the input should be
% specified as a cell array.
%
% See also FT_TIMELOCKGRANDAVERAGE, FT_FREQANALYSIS, FT_FREQDESCRIPTIVES

% FIXME averaging coherence is not possible if inputs contain different amounts of data (i.e. chan/freq/time)

% Copyright (C) 2005-2006, Robert Oostenveld
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

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'inputfile'),    cfg.inputfile  = [];         end
if ~isfield(cfg, 'outputfile'),   cfg.outputfile = [];         end

hasdata = nargin>1;
if ~isempty(cfg.inputfile) % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    for i=1:numel(cfg.inputfile)
      varargin{i} = loadvar(cfg.inputfile{i}, 'freq'); % read datasets from array inputfile
    end
  end
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq', 'feedback', 'no');
end

% set the defaults
if ~isfield(cfg, 'keepindividual'), cfg.keepindividual = 'no'; end
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';       end
if ~isfield(cfg, 'foilim'),         cfg.foilim = 'all';        end
if ~isfield(cfg, 'toilim'),         cfg.toilim = 'all';        end

% for backward compatibility with old data structures
if isfield(varargin{1}, 'sgn') && ~isfield(varargin{1}, 'label')
  warning('renaming "sng" field into label for backward compatibility');
  for s=1:Nsubj
    varargin{s}.label = varargin{s}.sgn;
  end
end

% for backward compatibility with old data structures
if isfield(varargin{1}, 'sgncmb') && ~isfield(varargin{1}, 'labelcmb')
  warning('renaming "sngcmb" field into labelcmb for backward compatibility');
  for s=1:Nsubj
    varargin{s}.labelcmb = varargin{s}.sgncmb;
  end
end

Nsubj    = length(varargin);
dimord   = varargin{1}.dimord;
haspow   = isfield(varargin{1}, 'powspctrm');             % this should always be true
hascoh   = isfield(varargin{1}, 'cohspctrm');
hasplv   = isfield(varargin{1}, 'plvspctrm');
hasfreq  = ~isempty(strfind(varargin{i}.dimord, 'freq')); % this should always be true
hastime  = ~isempty(strfind(varargin{i}.dimord, 'time'));
hasrpt   = ~isempty(strfind(varargin{i}.dimord, 'rpt'));
hastap   = ~isempty(strfind(varargin{i}.dimord, 'tap'));

% check whether the input data is suitable
if hasrpt
  error('the input data of each subject should be an average, use FREQDESCRIPTIVES first');
end
if hastap
  error('multiple tapers in the input are not supported');
end

if ischar(cfg.foilim) && strcmp(cfg.foilim, 'all')
  fbeg = -inf;
  fend =  inf;
else
  fbeg = cfg.foilim(1);
  fend = cfg.foilim(2);
end

if ischar(cfg.toilim) && strcmp(cfg.toilim, 'all')
  tbeg = -inf;
  tend =  inf;
else
  tbeg = cfg.toilim(1);
  tend = cfg.toilim(2);
end

% determine which channels, frequencies and latencies are available for all inputs
for i=1:Nsubj
  cfg.channel = ft_channelselection(cfg.channel, varargin{i}.label);
  if hasfreq
    fbeg = max(fbeg, varargin{i}.freq(1  ));
    fend = min(fend, varargin{i}.freq(end));
  end
  if hastime
    tbeg = max(tbeg, varargin{i}.time(1  ));
    tend = min(tend, varargin{i}.time(end));
  end
end
cfg.foilim = [fbeg fend];
cfg.toilim = [tbeg tend];

% select the data in all inputs
for i=1:Nsubj
  chansel = match_str(varargin{i}.label, cfg.channel);
  varargin{i}.label = varargin{i}.label(chansel);
  if hasfreq
    freqsel = nearest(varargin{i}.freq, fbeg):nearest(varargin{i}.freq, fend);
    varargin{i}.freq = varargin{i}.freq(freqsel);
  end
  if hastime
    timesel = nearest(varargin{i}.time, tbeg):nearest(varargin{i}.time, tend);
    varargin{i}.time = varargin{i}.time(timesel);
  end
  % select the overlapping samples in the power spectrum
  switch dimord
    case 'chan_freq'
      varargin{i}.powspctrm = varargin{i}.powspctrm(chansel,freqsel);
    case 'chan_freq_time'
      varargin{i}.powspctrm = varargin{i}.powspctrm(chansel,freqsel,timesel);
    case {'rpt_chan_freq' 'rpttap_chan_freq' 'subj_chan_freq'}
      varargin{i}.powspctrm = varargin{i}.powspctrm(:,chansel,freqsel);
    case {'rpt_chan_freq_time' 'rpttap_chan_freq_time' 'subj_chan_freq_time'}
      varargin{i}.powspctrm = varargin{i}.powspctrm(:,chansel,freqsel,timesel);
    otherwise
      error('unsupported dimord');
  end
end

% determine the size of the data to be averaged
dim = size(varargin{1}.powspctrm);
if hascoh,
  cohdim = size(varargin{1}.cohspctrm);
end
if hasplv,
  plvdim = size(varargin{1}.plvspctrm);
end

% give some feedback on the screen
if strcmp(cfg.keepindividual, 'no')
  if haspow, fprintf('computing average power over %d subjects\n', Nsubj); end
  if hascoh, fprintf('computing average coherence over %d subjects\n', Nsubj); end
  if hasplv, fprintf('computing average phase-locking value over %d subjects\n', Nsubj); end
else
  if haspow, fprintf('not computing grand average, but keeping individual power for %d subjects\n', Nsubj); end
  if hascoh, fprintf('not computing grand average, but keeping individual coherence for %d subjects\n', Nsubj); end
  if hasplv, fprintf('not computing grand average, but keeping individual phase-locking value for %d subjects\n', Nsubj); end
end

% allocate memory to hold the data
if strcmp(cfg.keepindividual, 'no')
  if haspow, s_pow = zeros(   dim); end
  if hascoh, s_coh = zeros(cohdim); end
  if hasplv, s_plv = zeros(plvdim); end
else
  if haspow, s_pow = zeros([Nsubj    dim]); end
  if hascoh, s_coh = zeros([Nsubj cohdim]); end
  if hasplv, s_plv = zeros([Nsubj plvdim]); end
end

for s=1:Nsubj
  if strcmp(cfg.keepindividual, 'no')
    % add this subject to the total sum
    if haspow, s_pow = s_pow + varargin{s}.powspctrm; end
    if hascoh, s_coh = s_coh + varargin{s}.cohspctrm; end
    if hasplv, s_plv = s_plv + varargin{s}.plvspctrm; end
  else
    % concatenate this subject to the rest
    if haspow, s_pow(s,:) = varargin{s}.powspctrm(:); end
    if hascoh, s_coh(s,:) = varargin{s}.cohspctrm(:); end
    if hasplv, s_plv(s,:) = varargin{s}.plvspctrm(:); end
  end
end

% collect the output data
grandavg.label = varargin{1}.label;
grandavg.freq  = varargin{1}.freq;
if isfield(varargin{1}, 'time')
  % remember the time axis
  grandavg.time = varargin{1}.time;
end
if isfield(varargin{1}, 'grad')
  warning('discarding gradiometer information because it cannot be averaged');
end
if isfield(varargin{1}, 'elec')
  warning('discarding electrode information because it cannot be averaged');
end

if strcmp(cfg.keepindividual, 'no')
  grandavg.dimord = dimord;
  grandavg.powspctrm = s_pow/Nsubj;
  if hascoh,
    grandavg.labelcmb = varargin{1}.labelcmb;
    grandavg.cohspctrm = s_coh/Nsubj;
  end
  if hasplv,
    grandavg.labelcmb = varargin{1}.labelcmb;
    grandavg.plvspctrm = s_plv/Nsubj;
  end
else
  grandavg.dimord = ['subj_' dimord];
  grandavg.powspctrm = s_pow;
  if hascoh,
    grandavg.labelcmb  = varargin{1}.labelcmb;
    grandavg.cohspctrm = s_coh;
  end
  if hasplv,
    grandavg.labelcmb  = varargin{1}.labelcmb;
    grandavg.plvspctrm = s_plv;
  end
end

cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.version.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();

% remember the configuration details of the input data
cfg.previous = [];
for i=1:length(varargin)
  try, cfg.previous{i} = varargin{i}.cfg; end
end

% remember the exact configuration details in the output
grandavg.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'freq', grandavg); % use the variable name "data" in the output file
end
