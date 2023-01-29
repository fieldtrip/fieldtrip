function [grandavg] = ft_freqgrandaverage(cfg, varargin)

% FT_FREQGRANDAVERAGE computes the average powerspectrum or time-frequency spectrum
% over multiple subjects
%
% Use as
%   [grandavg] = ft_freqgrandaverage(cfg, freq1, freq2, freq3...)
%
% The input data freq1..N are obtained from either FT_FREQANALYSIS with
% keeptrials=no or from FT_FREQDESCRIPTIVES. The configuration structure
% can contain
%   cfg.keepindividual = 'yes' or 'no' (default = 'no')
%   cfg.foilim         = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%   cfg.toilim         = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                        see FT_CHANNELSELECTION for details
%   cfg.parameter      = string or cell-array of strings indicating which
%                        parameter(s) to average. default is set to
%                        'powspctrm', if it is present in the data.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure. For this particular function, the input should be
% specified as a cell-array.
%
% See also FT_TIMELOCKGRANDAVERAGE, FT_FREQANALYSIS, FT_FREQDESCRIPTIVES,
% FT_FREQBASELINE

% FIXME averaging coherence is not possible if inputs contain different amounts of data (i.e. chan/freq/time)

% Copyright (C) 2005-2006, Robert Oostenveld
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
ft_preamble debug
ft_preamble loadvar varargin
ft_preamble provenance varargin

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq', 'feedback', 'no');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729

% set the defaults
cfg.keepindividual = ft_getopt(cfg, 'keepindividual', 'no');
cfg.channel        = ft_getopt(cfg, 'channel',    'all');
cfg.foilim         = ft_getopt(cfg, 'foilim',     'all');
cfg.toilim         = ft_getopt(cfg, 'toilim',     'all');
cfg.parameter      = ft_getopt(cfg, 'parameter',  []);

if isempty(cfg.parameter) && isfield(varargin{1}, 'powspctrm')
  cfg.parameter = 'powspctrm';
elseif isempty(cfg.parameter)
  ft_error('you should specify a valid parameter to average');
end

if ischar(cfg.parameter)
  cfg.parameter = {cfg.parameter};
end

Nsubj    = length(varargin);
dimord   = varargin{1}.dimord;
hasfreq  = ~isempty(strfind(varargin{i}.dimord, 'freq')); % this should always be true
hastime  = ~isempty(strfind(varargin{i}.dimord, 'time'));
hasrpt   = ~isempty(strfind(varargin{i}.dimord, 'rpt'));
hastap   = ~isempty(strfind(varargin{i}.dimord, 'tap'));

% check whether the input data is suitable
if hasrpt
  ft_error('the input data of each subject should be an average, use FT_FREQDESCRIPTIVES first');
end
if hastap
  ft_error('multiple tapers in the input are not supported');
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

% select the data in all inputs
for k=1:numel(cfg.parameter)
  
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
  
  % pick the selections
  for i=1:Nsubj
    if ~isfield(varargin{i}, cfg.parameter{k})
      ft_error('the field %s is not present in data structure %d', cfg.parameter{k}, i);
    end
    [dum, chansel] = match_str(cfg.channel, varargin{i}.label);
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
        varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(chansel,freqsel);
      case 'chan_freq_time'
        varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(chansel,freqsel,timesel);
      case {'rpt_chan_freq' 'rpttap_chan_freq' 'subj_chan_freq'}
        varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(:,chansel,freqsel);
      case {'rpt_chan_freq_time' 'rpttap_chan_freq_time' 'subj_chan_freq_time'}
        varargin{i}.(cfg.parameter{k}) = varargin{i}.(cfg.parameter{k})(:,chansel,freqsel,timesel);
      otherwise
        ft_error('unsupported dimord');
    end
  end % for i = subject
end % for k = parameter

% determine the size of the data to be averaged
dim = cell(1,numel(cfg.parameter));
for k=1:numel(cfg.parameter)
  dim{k} = size(varargin{1}.(cfg.parameter{k}));
end

% give some feedback on the screen
if strcmp(cfg.keepindividual, 'no')
  for k=1:numel(cfg.parameter)
    ft_info('computing average %s over %d subjects\n', cfg.parameter{k}, Nsubj);
  end
else
  for k=1:numel(cfg.parameter)
    ft_info('not computing average, but keeping individual %s for %d subjects\n', cfg.parameter{k}, Nsubj);
  end
end

% allocate memory to hold the data and collect it
for k=1:numel(cfg.parameter)
  if strcmp(cfg.keepindividual, 'no')
    tmp = zeros(dim{k});
    for s=1:Nsubj
      tmp = tmp + varargin{s}.(cfg.parameter{k})./Nsubj; % do a weighted running sum
    end
  elseif strcmp(cfg.keepindividual, 'yes')
    tmp = zeros([Nsubj dim{k}]);
    for s=1:Nsubj
      tmp(s,:,:,:,:) = varargin{s}.(cfg.parameter{k});
    end
  end
  grandavg.(cfg.parameter{k}) = tmp;
end

% collect the output data
grandavg.label = varargin{1}.label;
grandavg.freq  = varargin{1}.freq;
if isfield(varargin{1}, 'time')
  % remember the time axis
  grandavg.time = varargin{1}.time;
end
if isfield(varargin{1}, 'labelcmb')
  grandavg.labelcmb = varargin{1}.labelcmb;
end
if isfield(varargin{1}, 'grad')
  ft_warning('discarding gradiometer information because it cannot be averaged');
end
if isfield(varargin{1}, 'elec')
  ft_warning('discarding electrode information because it cannot be averaged');
end
if strcmp(cfg.keepindividual, 'yes')
  grandavg.dimord = ['subj_',varargin{1}.dimord];
elseif strcmp(cfg.keepindividual, 'no')
  grandavg.dimord = varargin{1}.dimord;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   varargin
ft_postamble provenance grandavg
ft_postamble history    grandavg
ft_postamble savevar    grandavg
