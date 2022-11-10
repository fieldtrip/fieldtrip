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
%   cfg.method         = string, 'across' or 'within' (default = 'across'), see below for details
%   cfg.parameter      = string, which parameter to average (default = 'avg')
%   cfg.channel        = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.latency        = [begin end] in seconds or 'all' (default = 'all')
%   cfg.keepindividual = string, 'yes' or 'no' (default = 'no')
%   cfg.nanmean        = string, can be 'yes' or 'no' (default = 'yes')
%   cfg.normalizevar   = string, 'N' or 'N-1' (default = 'N-1')
%
% If cfg.method = 'across', a plain average is performed, i.e. the requested
% parameter in each input argument is weighted equally in the average. This is useful
% when averaging across subjects. The variance-field will contain the variance across
% the parameter of interest, and the output dof-field will contain the number of
% input arguments.
%
% If cfg.method = 'within', a weighted average is performed, i.e. the requested
% parameter in each input argument is weighted according to the degrees of freedom in
% the dof-field. This is useful when averaging within subjects across blocks, e.g.
% when each block was recorded in a separate file. The variance-field will contain
% the variance across all input observations, and the output dof-field will contain
% the total number of observations.
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
    varargin{i} = rmfield(varargin{i}, 'trial');
    varargin{i}.dimord = 'chan_time';
    ft_warning('not using the trials, using the single-subject average to compute the grand average');
  else
    if isfield(varargin{i},'trial') && ~isfield(varargin{i},'avg')
      ft_error('input structure %d does not contain an average, use FT_TIMELOCKANALYSIS first', i);
    end
  end
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'no');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729

% set the defaults
cfg.method         = ft_getopt(cfg, 'method'        , 'across');
cfg.parameter      = ft_getopt(cfg, 'parameter'     , 'avg');
cfg.channel        = ft_getopt(cfg, 'channel'       , 'all');
cfg.latency        = ft_getopt(cfg, 'latency'       , 'all');
cfg.keepindividual = ft_getopt(cfg, 'keepindividual', 'no');
cfg.nanmean        = ft_getopt(cfg, 'nanmean'       , 'yes');
cfg.normalizevar   = ft_getopt(cfg, 'normalizevar'  , 'N-1');

if iscell(cfg.parameter)
  if numel(cfg.parameter)>1
    ft_error('only a single parameter can be specified to be averaged');
  else
    cfg.parameter = cfg.parameter{1};
  end
end

if strcmp(cfg.parameter, 'trial')
  ft_error('not supporting averaging over the repetition dimension, please use FT_TIMELOCKANALYSIS');
end

% select trials and channels of interest
orgcfg = cfg;
tmpcfg = keepfields(cfg, {'parameter', 'channel', 'tolerance', 'latency', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});
% do not use the default option returned by FT_SELECTDATA, but the original one for this function
cfg.nanmean = orgcfg.nanmean;

% determine the size of the data to be averaged
datsiz = size(varargin{1}.(cfg.parameter));
dimord = getdimord(varargin{1}, cfg.parameter);
nsubj  = length(varargin);

% whether to normalize the variance with N or N-1, see VAR
normalizewithN = strcmpi(cfg.normalizevar, 'N');

% whether nans should persist in the output or be treated as missing values
if istrue(cfg.nanmean)
  mymean = @nanmean;
  myvar  = @nanvar;
  mysum  = @nansum;
else
  mymean = @mean;
  myvar  = @var;
  mysum  = @sum;
end

if istrue(cfg.keepindividual)
  fprintf('not computing average, but keeping individual %s for %d subjects\n', cfg.parameter, nsubj);
  
  % allocate memory to hold the data and collect it
  dat = zeros([nsubj, datsiz]);
  for s=1:nsubj
    dat(s, :, :, :) = varargin{s}.(cfg.parameter);
  end
  grandavg.individual = dat; % Nsubj x Nchan x Nsamples
  
else
  dat = nan([nsubj, datsiz]);
  dof = nan([nsubj, datsiz]);
  
  switch cfg.method
    case 'across'
      fprintf('computing average of %s across %d subjects\n', cfg.parameter, nsubj);
      for s=1:nsubj
        dat(s, :, :, :) = varargin{s}.(cfg.parameter);
      end
      
      % compute the mean and variance across subjects
      grandavg.avg = reshape(mymean(dat, 1),                 datsiz);
      grandavg.var = reshape(myvar(dat, normalizewithN, 1),  datsiz);
      grandavg.dof = reshape(sum(isfinite(dat), 1),          datsiz);
      
      if normalizewithN
        % just to be sure
        grandavg.var(grandavg.dof<=0) = NaN;
      else
        % see https://stats.stackexchange.com/questions/4068/how-should-one-define-the-sample-variance-for-scalar-input
        % the fieldtrip/external/stats/nanvar implementation behaves differently here than Mathworks VAR and NANVAR implementations
        grandavg.var(grandavg.dof<=1) = NaN;
      end
      
    case 'within'
      fprintf('computing average of %s within subjects over %d blocks\n', cfg.parameter, nsubj);
      
      for s=1:nsubj
        dat(s, :, :, :) = varargin{s}.(cfg.parameter);
        dof(s, :, :, :) = varargin{s}.dof;
      end % for nsub
      
      % compute the weighted mean across input arguments
      grandavg.avg = mysum(dat .* dof, 1) ./ sum(dof, 1);
      grandavg.avg = reshape(grandavg.avg, datsiz);
      grandavg.dof = sum(dof, 1);
      grandavg.dof = reshape(grandavg.dof, datsiz);
      
      % this follows https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance with frequency weights
      if normalizewithN
        grandavg.var = mysum(dof .* (dat - repmat(reshape(grandavg.avg, [1 datsiz]), [nsubj 1 1])).^2, 1) ./  sum(dof, 1);
        grandavg.var = reshape(grandavg.var, datsiz);
      else
        grandavg.var = mysum(dof .* (dat - repmat(reshape(grandavg.avg, [1 datsiz]), [nsubj 1 1])).^2, 1) ./ (sum(dof, 1) - 1);
        grandavg.var = reshape(grandavg.var, datsiz);
      end
      
    otherwise
      ft_error('unsupported value for cfg.method')
      
  end % switch
end % if keepindividual

% collect the output data
grandavg = copyfields(varargin{1}, grandavg, {'time', 'freq', 'label', 'labelcmb'});

% sensor positions should only be present in the output when they are identical in all inputs
hasgrad = cellfun(@(x) isfield(x, 'grad'), varargin(:));
if all(hasgrad) % check if positions are different between subjects
  samegrad = cellfun(@(x) isequal(varargin{1}.grad, x.grad), varargin(2:end));
  if all(samegrad)
    grandavg.grad = varargin{1}.grad;
  else
    ft_warning('discarding gradiometer information because it is not identical in all inputs');
  end
end

haselec = cellfun(@(x) isfield(x, 'elec'), varargin(:));
if all(haselec) % check if positions are different between subjects
  sameelec = cellfun(@(x) isequal(varargin{1}.elec, x.elec), varargin(2:end));
  if all(sameelec)
    grandavg.elec = varargin{1}.elec;
  else
    ft_warning('discarding electrode information because it is not identical in all inputs');
  end
end

hasopto = cellfun(@(x) isfield(x, 'opto'), varargin(:));
if all(hasopto) % check if positions are different between subjects
  sameopto = cellfun(@(x) isequal(varargin{1}.opto, x.opto), varargin(2:end));
  if all(sameopto)
    grandavg.opto = varargin{1}.opto;
  else
    ft_warning('discarding optode information because it is not identical in all inputs');
  end
end

% set dimord
if strcmp(cfg.keepindividual, 'yes')
  grandavg.dimord = ['subj_', dimord];
elseif strcmp(cfg.keepindividual, 'no')
  grandavg.dimord =  dimord;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   varargin
ft_postamble provenance grandavg
ft_postamble history    grandavg
ft_postamble savevar    grandavg
