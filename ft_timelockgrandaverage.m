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
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'),
%                        see FT_CHANNELSELECTION for details
%   cfg.latency        = [begin end] in seconds or 'all' (default = 'all')
%   cfg.keepindividual = 'yes' or 'no' (default = 'no')
%   cfg.normalizevar   = 'N' or 'N-1' (default = 'N-1')
%   cfg.method         = 'across' (default) or 'within', see below.
%   cfg.parameter      = string or cell-array indicating which
%                        parameter to average. default is set to
%                        'avg', if it is present in the data.
%
% If cfg.method = 'across', a plain average is performed, i.e. the
% requested parameter in each input argument is weighted equally in the
% average. This is useful when averaging across subjects. The
% variance-field will contain the variance across the parameter of
% interest, and the dof-field will contain the number of input arguments.
%
% If cfg.method = 'within', a weighted average is performed, i.e. the
% requested parameter in each input argument is weighted according to the
% dof-field. This is useful when averaging across blocks within subjects.
% The variance-field will contain the variance across all input
% observations, and the dof-field will contain the number of observations.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure. For this particular function, the input should be
% structured as a cell-array.
%
% See also FT_TIMELOCKANALYSIS, FT_TIMELOCKSTATISTICS, FT_TIMELOCKBASELINE

% Copyright (C) 2003-2006, Jens Schwarzbach
% Copyright (C) 2013, Burkhard Maess
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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% return immediately after distributed execution
if ~isempty(ft_getopt(cfg, 'distribute'))
  return
end

% check if the input data is valid for this function
for i=1:length(varargin)
  if isfield(varargin{i},'trial') && isfield(varargin{i},'avg') % see bug2372 (dieloz)
    varargin{i} = rmfield(varargin{i},'trial');
    ft_warning('depreciating trial field: using the avg to compute the grand average');
    if strcmp(varargin{i}.dimord,'rpt_chan_time')
      varargin{i}.dimord = 'chan_time';
    end
  else
    if isfield(varargin{i},'trial') && ~isfield(varargin{i},'avg')
      ft_error('input dataset %d does not contain avg field: see ft_timelockanalysis', i);
    end
  end
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

if strcmp(cfg.parameter,'trial')
  ft_error('not supporting averaging over the repetition dimension');
end

Nsubj    = length(varargin);
dimord   = getdimord(varargin{1}, cfg.parameter);
hastime  = contains(dimord, 'time');
hasdof   = isfield(varargin{1}, 'dof');

if ischar(cfg.latency) && strcmp(cfg.latency, 'all')
  tbeg = min(varargin{1}.time);
  tend = max(varargin{1}.time);
else
  tbeg = cfg.latency(1);
  tend = cfg.latency(2);
end

% ensure that the data in all inputs has the same channels, time-axis, etc.
tmpcfg = keepfields(cfg, {'parameter', 'channel', 'latency', 'showcallinfo'});
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

% determine the size of the data to be averaged
datsiz = size(varargin{1}.(cfg.parameter));

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
avgmat = zeros([Nsubj, datsiz]);

if strcmp(cfg.keepindividual, 'yes')
  for s=1:Nsubj
    avgmat(s, :, :) = varargin{s}.(cfg.parameter);
  end
  grandavg.individual = avgmat; % Nsubj x Nchan x Nsamples
  
else % ~strcmp(cfg.keepindividual, 'yes')
  avgdof  = ones([Nsubj, datsiz]);
  avgvar  = zeros([Nsubj, datsiz]);
  for s=1:Nsubj
    switch cfg.method
      case 'across'
        avgmat(s, :, :, :) = varargin{s}.(cfg.parameter);
        avgvar(s, :, :, :) = varargin{s}.(cfg.parameter) .^2; % preparing the computation of the variance
      case 'within'
        avgmat(s, :, :, :) = varargin{s}.(cfg.parameter).*varargin{s}.dof;
        avgdof(s, :, :, :) = varargin{s}.dof;
        if min(varargin{s}.dof(:))>1 % otherwise the variance is not valid
          avgvar(s, :, :, :) = (varargin{s}.(cfg.parameter).^2).*varargin{s}.dof;
          % avgvar(s, :, :, :) = varargin{s}.var .* (varargin{s}.dof-1); % reversing the last div in ft_timelockanalysis
        else
          avgvar(s, :, :, :) = zeros(datsiz); % shall we remove the .var field from the structure under these conditions ?
        end
      otherwise
        ft_error('unsupported value for cfg.method')
    end % switch
  end
  % average across subject dimension
  ResultDOF      = reshape(sum(avgdof, 1), datsiz);
  grandavg.avg   = reshape(sum(avgmat, 1), datsiz)./ResultDOF; % computes both means (plain and weighted)
  % Nchan x Nsamples, skips the singleton
  % if strcmp(cfg.method, 'across')
  ResultVar      = reshape(sum(avgvar,1), datsiz)-reshape(sum(avgmat,1), datsiz).^2./ResultDOF;
  % else  % cfg.method = 'within'
  % ResultVar      = reshape(sum(avgvar, 1), datsiz); % subtraction of means was done for each block already
  % end
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

switch cfg.method
  
  case 'across'
    if isfield(varargin{1}, 'grad') % positions are different between subjects
      ft_warning('discarding gradiometer information because it cannot be averaged');
    end
    if isfield(varargin{1}, 'elec') % positions are different between subjects
      ft_warning('discarding electrode information because it cannot be averaged');
    end
    
  case 'within'
    % misses the test for all equal grad fields (should be the case for
    % averaging across blocks, if all block data is corrected for head
    % position changes
    if isfield(varargin{1}, 'grad')
      grandavg.grad = varargin{1}.grad;
    end
    if isfield(varargin{1}, 'elec')
      grandavg.elec = varargin{1}.elec;
    end
    
  otherwise
    ft_error('unsupported method "%s"', cfg.method);
end

if strcmp(cfg.keepindividual, 'yes')
  grandavg.dimord = ['subj_', dimord];
else
  grandavg.dimord = dimord;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance grandavg
ft_postamble history    grandavg
ft_postamble savevar    grandavg
