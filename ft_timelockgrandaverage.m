function [grandavg] = ft_timelockgrandaverage(cfg, varargin)

% FT_TIMELOCKGRANDAVERAGE computes ERF/ERP average and variance
% over multiple subjects or over blocks within one subject
%
% Use as
%   [grandavg] = ft_timelockgrandaverage(cfg, avg1, avg2, avg3, ...)
%
% where
%   avg1..N are the ERF/ERP averages as obtained from FT_TIMELOCKANALYSIS
%
% and cfg is a configuration structure with
%  cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                       see FT_CHANNELSELECTION for details
%  cfg.latency        = [begin end] in seconds or 'all' (default = 'all')
%  cfg.keepindividual = 'yes' or 'no' (default = 'no')
%  cfg.normalizevar   = 'N' or 'N-1' (default = 'N-1')
%  cfg.method         = 'across' (default) or 'within', see below.
%  cfg.parameter      = string or cell-array indicating which
%                        parameter to average. default is set to
%                        'avg', if it is present in the data.
%
% If cfg.method = 'across', an plain average is performed, i.e. the
% requested parameter in each input argument is weighted equally in the
% average. This is useful when averaging across subjects. The
% variance-field will contain the variance across the parameter of
% interest, and the dof-field will contain the number of input arguments.
% If cfg.method = 'within', a weighted average is performed, i.e. the
% requested parameter in each input argument is weighted according to the
% dof-field. This is useful when averaging across blocks within subjects.
% The variance-field will contain the variance across all input
% observations, and the dof-field will contain the number of observations.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure. For this particular function, the input should be
% structured as a cell array.
%
% See also FT_TIMELOCKANALYSIS, FT_TIMELOCKSTATISTICS

% Copyright (C) 2003-2006, Jens Schwarzbach
% Copyright (C) 2013, Burkhard Maess
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
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar varargin

% return immediately after distributed execution
if ~isempty(ft_getopt(cfg, 'distribute'))
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'no');
end

% set the defaults
cfg.keepindividual = ft_getopt(cfg, 'keepindividual', 'no');
cfg.channel        = ft_getopt(cfg, 'channel',    'all');
cfg.latency        = ft_getopt(cfg, 'latency',    'all');
cfg.normalizevar   = ft_getopt(cfg, 'normalizevar', 'N-1');
cfg.method         = ft_getopt(cfg, 'method',    'across');
cfg.parameter      = ft_getopt(cfg, 'parameter',  'avg');

if iscell(cfg.parameter)
  if numel(cfg.parameter)>1
    fprintf('ft_grandaverage supports a single parameter to be averaged, taking the first parameter from the specified list: %s\n', cfg.parameter{1});
  end
  cfg.parameter = cfg.parameter{1};
end

Nsubj    = length(varargin);
dimord   = varargin{1}.dimord;
hastime  = ~isempty(strfind(varargin{1}.dimord, 'time'));
hasrpt   = ~isempty(strfind(varargin{1}.dimord, 'rpt'));

% check whether the input data is suitable
if hasrpt
  fprintf('ignoring the single-subject repetition dimension\n');
  if strcmp(cfg.parameter, 'trial')
    error('not supporting averaging over the repetition dimension');
  end
end

if ischar(cfg.latency) && strcmp(cfg.latency, 'all')
  tbeg = min(varargin{1}.time);
  tend = max(varargin{1}.time);
else
  tbeg = cfg.latency(1);
  tend = cfg.latency(2);
end

% select the data in all inputs
% determine which channels, and latencies are available for all inputs
for i=1:Nsubj
  cfg.channel = ft_channelselection(cfg.channel, varargin{i}.label);
  if hastime
    tbeg = max(tbeg, varargin{i}.time(1  ));
    tend = min(tend, varargin{i}.time(end));
  end
end
cfg.latency = [tbeg tend];

% pick the selections
for i=1:Nsubj
  if ~isfield(varargin{i}, cfg.parameter)
    error('the field %s is not present in data structure %d', cfg.parameter, i);
  end
  chansel = match_str(varargin{i}.label, cfg.channel);
  varargin{i}.label = varargin{i}.label(chansel);
  if hastime
    timesel = nearest(varargin{i}.time, cfg.latency(1)):nearest(varargin{i}.time, cfg.latency(2));
    varargin{i}.time = varargin{i}.time(timesel);
  end
  % select the overlapping samples in the data matrix
  switch dimord
    case 'chan'
      varargin{i}.(cfg.parameter) = varargin{i}.(cfg.parameter)(chansel);
    case 'chan_time'
      varargin{i}.(cfg.parameter) = varargin{i}.(cfg.parameter)(chansel,timesel);
    case {'rpt_chan' 'subj_chan'}
      varargin{i}.(cfg.parameter) = varargin{i}.(cfg.parameter)(chansel);
      varargin{i}.dimord = 'chan';
    case {'rpt_chan_time' 'subj_chan_time'}
      varargin{i}.(cfg.parameter) = varargin{i}.(cfg.parameter)(chansel,timesel);
      varargin{i}.dimord = 'chan_time';
    otherwise
      error('unsupported dimord');
  end
end % for i = subject

% determine the size of the data to be averaged
%dim = cell(1,numel(cfg.parameter));
dim{1} = size(varargin{1}.(cfg.parameter));

% give some feedback on the screen
if strcmp(cfg.keepindividual, 'no')
  if strcmp(cfg.method, 'across')
    fprintf('computing average of %s over %d subjects\n', cfg.parameter, Nsubj);
  else % cfg.method = 'within'
    fprintf('computing average of %s over %d blocks\n', cfg.parameter, Nsubj);
  end
else
  fprintf('not computing average, but keeping individual %s for %d subjects\n', cfg.parameter, Nsubj);
end

% allocate memory to hold the data and collect it
avgmat = zeros([Nsubj, dim{1}]);
if strcmp(cfg.keepindividual, 'yes')
  for s=1:Nsubj
    avgmat(s, :, :) = varargin{s}.(cfg.parameter);
  end
  grandavg.individual = avgmat;         % Nsubj x Nchan x Nsamples
else % ~strcmp(cfg.keepindividual, 'yes')
  avgdof  = ones([Nsubj, dim{1}]);
  avgvar  = zeros([Nsubj, dim{1}]);
  for s=1:Nsubj
    if strcmp(cfg.method, 'across')
      avgmat(s, :, :) = varargin{s}.(cfg.parameter);
      avgvar(s, :, :) = varargin{s}.(cfg.parameter) .^2;     % preparing the computation of the variance
    else % cfg.method = 'within'
      avgmat(s, :, :) = varargin{s}.(cfg.parameter).*varargin{s}.dof;
      avgdof(s, :, :) = varargin{s}.dof;
      if min(varargin{s}.dof(:))>1 % otherwise the variance is not valid
        avgvar(s, :, :) = (varargin{s}.(cfg.parameter).^2).*varargin{s}.dof;
        % avgvar(s, :, :) = varargin{s}.var .* (varargin{s}.dof-1); % reversing the last div in ft_timelockanalysis
      else
        avgvar(s, :, :) = zeros([dim{1}]); % shall we remove the .var field from the structure under these conditions ?
      end
    end
  end
  % average across subject dimension
  ResultDOF      = reshape(sum(avgdof, 1), dim{1});
  grandavg.avg   = reshape(sum(avgmat, 1), dim{1})./ResultDOF; % computes both means (plain and weighted)
  % Nchan x Nsamples, skips the singleton
  %if strcmp(cfg.method, 'across')
  ResultVar      = reshape(sum(avgvar,1), dim{1})-reshape(sum(avgmat,1), dim{1}).^2./ResultDOF;
  %else  % cfg.method = 'within'
    %ResultVar      = squeeze(sum(avgvar,1)); % subtraction of means was done for each block already
  %end
  switch cfg.normalizevar
    case 'N-1'
      ResultVar = ResultVar./(ResultDOF-1);
    case 'N'
      ResultVar = ResultVar./ResultDOF;
  end
  grandavg.var = ResultVar;
  grandavg.dof = ResultDOF;
end

% collect the output data
if isfield(varargin{1}, 'time')
  % remember the time axis
  grandavg.time = varargin{1}.time;
end
grandavg.label = varargin{1}.label;

if isfield(varargin{1}, 'labelcmb')
  grandavg.labelcmb = varargin{1}.labelcmb;
end

if strcmp(cfg.method, 'across')
  if isfield(varargin{1}, 'grad') % positions are different between subjects
    warning('discarding gradiometer information because it cannot be averaged');
  end
  if isfield(varargin{1}, 'elec') % positions are different between subjects
    warning('discarding electrode information because it cannot be averaged');
  end
else  % cfg.method = 'within'
  % misses the test for all equal grad fields (should be the case for
  % averaging across blocks, if all block data is corrected for head
  % position changes
  if isfield(varargin{1}, 'grad')
    grandavg.grad = varargin{1}.grad;
  end
  if isfield(varargin{1}, 'elec')
    grandavg.elec = varargin{1}.elec;
  end
  
end
if strcmp(cfg.keepindividual, 'yes')
  grandavg.dimord = ['subj_',varargin{1}.dimord];
elseif strcmp(cfg.keepindividual, 'no')
  grandavg.dimord = varargin{1}.dimord;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous varargin
ft_postamble history grandavg
ft_postamble savevar grandavg
