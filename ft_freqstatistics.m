function [stat] = ft_freqstatistics(cfg, varargin)

% FT_FREQSTATISTICS computes significance probabilities and/or critical
% values of a parametric statistical test or a non-parametric permutation
% test.
%
% Use as
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
% where the input data is the result from FT_FREQANALYSIS, FT_FREQDESCRIPTIVES
% or from FT_FREQGRANDAVERAGE.
%
% The configuration can contain the following options for data selection
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
%   cfg.latency     = [begin end] in seconds or 'all' (default = 'all')
%   cfg.trials      = trials to be included or 'all'  (default = 'all')
%   cfg.frequency   = [begin end], can be 'all'       (default = 'all')
%   cfg.avgoverchan = 'yes' or 'no'                   (default = 'no')
%   cfg.avgovertime = 'yes' or 'no'                   (default = 'no')
%   cfg.avgoverfreq = 'yes' or 'no'                   (default = 'no')
%   cfg.parameter   = string                          (default = 'powspctrm')
%
% If you specify cfg.correctm='cluster', then the following is required
%   cfg.neighbours  = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
%
% Furthermore, the configuration should contain
%   cfg.method       = different methods for calculating the significance probability and/or critical value
%                    'montecarlo'    get Monte-Carlo estimates of the significance probabilities and/or critical values from the permutation distribution,
%                    'analytic'      get significance probabilities and/or critical values from the analytic reference distribution (typically, the sampling distribution under the null hypothesis),
%                    'stats'         use a parametric test from the Matlab statistics toolbox,
%                    'crossvalidate' use crossvalidation to compute predictive performance
%
%   cfg.design       = Nxnumobservations: design matrix (for examples/advice, please see the Fieldtrip wiki,
%                      especially cluster-permutation tutorial and the 'walkthrough' design-matrix section)
%
% The other cfg options depend on the method that you select. You
% should read the help of the respective subfunction FT_STATISTICS_XXX
% for the corresponding configuration options and for a detailed
% explanation of each method.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_FREQANALYSIS, FT_FREQDESCRIPTIVES, FT_FREQGRANDAVERAGE

% TODO change cfg.frequency in all functions to cfg.foi or cfg.foilim

% Copyright (C) 2005-2010, Robert Oostenveld
% Copyright (C) 2011, Jan-Mathijs Schoffelen
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar varargin

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',     {'approach',   'method'});
cfg = ft_checkconfig(cfg, 'required',    {'method'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'transform'});

% set the defaults
cfg.outputfile  = ft_getopt(cfg, 'outputfile',  []);
cfg.parameter   = ft_getopt(cfg, 'parameter',   []); % the default is assigned further down
cfg.channel     = ft_getopt(cfg, 'channel',     'all');
cfg.latency     = ft_getopt(cfg, 'latency',     'all');
cfg.trials      = ft_getopt(cfg, 'trials',      'all');
cfg.frequency   = ft_getopt(cfg, 'frequency',   'all');
cfg.avgoverchan = ft_getopt(cfg, 'avgoverchan', 'no');
cfg.avgoverfreq = ft_getopt(cfg, 'avgoverfreq', 'no');
cfg.avgovertime = ft_getopt(cfg, 'avgovertime', 'no');
cfg.design      = ft_getopt(cfg, 'design',      '');

% get the design from the information in cfg and data.
if ~isfield(cfg,'design') || isempty(cfg.design)
  error('you should provide a design matrix in the cfg');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data bookkeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndata = numel(varargin);

% set the parameter to default powspctrm only if present in the data
if isempty(cfg.parameter) && isfield(varargin{1}, 'powspctrm')
  cfg.parameter = 'powspctrm';
elseif isempty(cfg.parameter)
  error('You need to specify a cfg.parameter, because the default (powspctrm) is not present in the input data');
end

% check if the input data is valid for this function
hastime  = false(Ndata,1);
hasparam = false(Ndata,1);
for i=1:Ndata
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq', 'feedback', 'no');
  hastime(i)  = isfield(varargin{i}, 'time');
  haschan(i)  = ~isempty(strfind(varargin{i}.dimord, 'chan')) ...
                && numel(varargin{i}.label) > 1;
  hasparam(i) = isfield(varargin{i}, cfg.parameter);
end

if sum(hastime)~=Ndata && sum(hastime)~=0
  error('the input data structures should either all contain a time axis, or none of them should');
end
hastime = sum(hastime)==Ndata;

if sum(hasparam)~=Ndata
  error('the input data structures should all contain the parameter %s', cfg.parameter);
end

haschan = sum(haschan)==Ndata;

% check whether channel neighbourhood information is needed and present
if isfield(cfg, 'correctm') && strcmp(cfg.correctm, 'cluster') ...
    && ~strcmp(cfg.avgoverchan, 'yes') && haschan
  cfg = ft_checkconfig(cfg, 'required', {'neighbours'});    
end

% get frequency, latency and channels which are present in all subjects
fmin = -inf;
fmax =  inf;
tmin = -inf;
tmax =  inf;
for i=1:Ndata
  fmin = max(fmin, varargin{i}.freq(1));
  fmax = min(fmax, varargin{i}.freq(end));
  if hastime
    tmin = max(tmin, varargin{i}.time(1));
    tmax = min(tmax, varargin{i}.time(end));
  end
  if i==1
    % FIXME deal with channelcmb
    if isfield(varargin{i}, 'labelcmb')
      error('support for data containing linearly indexed bivariate quantities, i.e. containing a ''labelcmb'' is not yet implemented');
    else
      chan = varargin{i}.label;
    end
  else
    chan = varargin{i}.label(ismember(varargin{i}.label, chan));
  end
end

if ischar(cfg.frequency) && strcmp(cfg.frequency, 'all')
  cfg.frequency = [fmin fmax];
elseif ischar(cfg.frequency)
  error('unsupported value for ''cfg.frequency''');
end

% overrule user-specified settings
cfg.frequency = [max(cfg.frequency(1), fmin), min(cfg.frequency(2), fmax)];
fprintf('computing statistic over the frequency range [%1.3f %1.3f]\n', cfg.frequency(1), cfg.frequency(2));

if hastime
  if ischar(cfg.latency) && strcmp(cfg.latency, 'all')
    cfg.latency = [tmin tmax];
  elseif ischar(cfg.latency)
    error('unsupported value for ''cfg.latency''');
  end
  
  % overrule user-specified settings
  cfg.latency = [max(cfg.latency(1), tmin), min(cfg.latency(2), tmax)];
  fprintf('computing statistic over the time range [%1.3f %1.3f]\n', cfg.latency(1), cfg.latency(2));
end

% only do those channels present in the data
cfg.channel = ft_channelselection(cfg.channel, chan);

if ~ischar(cfg.trials)
  if Ndata==1
    varargin{1} = ft_selectdata(varargin{1}, 'rpt', cfg.trials);
  else
    error('subselection of trials is only allowed with a single data structure as input');
  end
end

% intersect the data and combine it into one structure
if hastime
  if haschan
    data =  ft_selectdata(varargin{:}, 'param', cfg.parameter, 'avgoverrpt', false, ...
      'toilim', cfg.latency, 'avgovertime', cfg.avgovertime, ...
      'foilim',  cfg.frequency, 'avgoverfreq', cfg.avgoverfreq, ...
      'channel', cfg.channel, 'avgoverchan', cfg.avgoverchan);
  else
    data =  ft_selectdata(varargin{:}, 'param', cfg.parameter, 'avgoverrpt', false, ...
    'toilim', cfg.latency, 'avgovertime', cfg.avgovertime, ...
    'foilim',  cfg.frequency, 'avgoverfreq', cfg.avgoverfreq, ...
    'avgoverchan', cfg.avgoverchan);
  end
else
  if haschan
    data =  ft_selectdata(varargin{:}, 'param', cfg.parameter, 'avgoverrpt', false, ...
      'foilim',  cfg.frequency, 'avgoverfreq', cfg.avgoverfreq, ...
      'channel', cfg.channel, 'avgoverchan', cfg.avgoverchan);
  else
    data =  ft_selectdata(varargin{:}, 'param', cfg.parameter, 'avgoverrpt', false, ...
      'foilim',  cfg.frequency, 'avgoverfreq', cfg.avgoverfreq, ...
      'avgoverchan', cfg.avgoverchan);
  end
end

% keep the sensor info, just in case
if isfield(varargin{1}, 'elec')
  data.elec = varargin{1}.elec;
elseif isfield(varargin{1}, 'grad')
  data.grad = varargin{1}.grad;
end

% ensure that we don't touch this any more
clear varargin;

% create the 'dat' matrix here
dat        = data.(cfg.parameter);
siz        = size(dat);
dimtok     = tokenize(data.dimord, '_');
% check for occurence of channel dimension
chandim     = find(ismember(dimtok, {'chan'}));
if isempty(chandim)
  dimtok(3:end+1) = dimtok(2:end);
  dimtok{2} = 'chan';
  siz = [siz(1) 1 siz(2:end)];
  dat = reshape(dat, siz);
end
rptdim     = find(ismember(dimtok, {'rpt' 'subj' 'rpttap'}));
permutevec = [setdiff(1:numel(siz), rptdim) rptdim];       % permutation vector to put the repetition dimension as last dimension
reshapevec = [prod(siz(permutevec(1:end-1))) siz(rptdim) 1]; % reshape vector to reshape into 2D
dat        = reshape(permute(dat, permutevec), reshapevec);% actually reshape the data

reduceddim = setdiff(1:numel(siz), rptdim);
cfg.dim    = [siz(reduceddim) 1];   % store dimensions of the output of the statistics function in the cfg
cfg.dimord = '';
for k = 1:numel(reduceddim)
  cfg.dimord = [cfg.dimord, '_', dimtok{reduceddim(k)}];
end
cfg.dimord = cfg.dimord(2:end); % store the dimord of the output in the cfg

if size(cfg.design,2)~=size(dat,2)
  error('the number of observations in the design does not match the number of observations in the data');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the statistic, using the data-independent statistical subfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine the function handle to the intermediate-level statistics function
if exist(['ft_statistics_' cfg.method], 'file')
  statmethod = str2func(['ft_statistics_' cfg.method]);
else
  error('could not find the corresponding function for cfg.method="%s"\n', cfg.method);
end
fprintf('using "%s" for the statistical testing\n', func2str(statmethod));

% check that the design completely describes the data
if size(dat,2) ~= size(cfg.design,2)
  error('the size of the design matrix does not match the number of observations in the data');
end

% determine the number of output arguments
num = nargout(statmethod);

% perform the statistical test
if num>1
  [stat, cfg] = statmethod(cfg, dat, cfg.design);
else
  [stat] = statmethod(cfg, dat, cfg.design);
end

if isstruct(stat)
  % the statistical output contains multiple elements, e.g. F-value, beta-weights and probability
  statfield = fieldnames(stat);
else
  % only the probability was returned as a single matrix, reformat into a structure
  dum = stat; stat = []; % this prevents a Matlab warning that appears from release 7.0.4 onwards
  stat.prob = dum;
  statfield = fieldnames(stat);
end

% post-process the output data
haschan     = isfield(data, 'label');    % this one remains relevant, even after averaging over channels
haschancmb  = isfield(data, 'labelcmb'); % this one remains relevant, even after averaging over channels
hasfreq     = strcmp(cfg.avgoverfreq, 'no') && ~any(isnan(data.freq));
hastime     = hastime && strcmp(cfg.avgovertime, 'no');
stat.dimord = '';

if haschan
  stat.label  = data.label;
elseif haschancmb
  stat.labelcmb = data.labelcmb;
end

if hasfreq
  stat.freq = data.freq;
end

if hastime
  stat.time = data.time;
end

if ~isempty(stat.dimord)
  % remove the last '_'
  stat.dimord = stat.dimord(1:(end-1));
end

for i=1:length(statfield)
  if numel(stat.(statfield{i})) == prod(cfg.dim)
    % reshape the fields that have the same dimension as the input data
    stat.(statfield{i}) = reshape(stat.(statfield{i}), cfg.dim);
  end
end

% FIXME squeeze out the appropriate dimords if avgoverfreq etc.
stat.dimord = cfg.dimord;

% HACK if a bivariate statistic is in the output, replace label with the
% appropriate labelcmb
if isfield(cfg,'statistic') && strcmp(cfg.statistic, 'indepsamplesZcoh') && isfield(stat, 'label')
  stat.labelcmb(:,1) = stat.label(cfg.chancmbindx(:,1));
  stat.labelcmb(:,2) = stat.label(cfg.chancmbindx(:,2));
  stat = rmfield(stat, 'label');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous varargin
ft_postamble history stat
ft_postamble savevar stat



