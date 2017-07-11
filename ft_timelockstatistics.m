function [stat] = ft_timelockstatistics(cfg, varargin)

% FT_TIMELOCKSTATISTICS  computes significance probabilities and/or critical values of a parametric statistical test
% or a non-parametric permutation test.
%
% Use as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
% where the input data is the result from either FT_TIMELOCKANALYSIS or
% FT_TIMELOCKGRANDAVERAGE.
%
% The configuration can contain the following options for data selection
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
%   cfg.latency     = [begin end] in seconds or 'all' (default = 'all')
%   cfg.avgoverchan = 'yes' or 'no'                   (default = 'no')
%   cfg.avgovertime = 'yes' or 'no'                   (default = 'no')
%   cfg.parameter   = string                          (default = 'trial' or 'avg')
%
% Furthermore, the configuration should contain
%   cfg.method       = different methods for calculating the significance probability and/or critical value
%                    'montecarlo'    get Monte-Carlo estimates of the significance probabilities and/or critical values from the permutation distribution,
%                    'analytic'      get significance probabilities and/or critical values from the analytic reference distribution (typically, the sampling distribution under the null hypothesis),
%                    'stats'         use a parametric test from the MATLAB statistics toolbox,
%                    'crossvalidate' use crossvalidation to compute predictive performance
%
% The other cfg options depend on the method that you select. You
% should read the help of the respective subfunction FT_STATISTICS_XXX
% for the corresponding configuration options and for a detailed
% explanation of each method.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_TIMELOCKANALYSIS, FT_TIMELOCKGRANDAVERAGE

% Copyright (C) 2005-2014, Robert Oostenveld
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required',    {'method', 'design'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'trials'}); % this used to be present until 24 Dec 2014, but was deemed too confusing by Robert

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'timelock', 'feedback', 'no');
end

% set the defaults
cfg.parameter   = ft_getopt(cfg, 'parameter'); % default is set below
cfg.correctm    = ft_getopt(cfg, 'correctm');
cfg.channel     = ft_getopt(cfg, 'channel',     'all');
cfg.avgoverchan = ft_getopt(cfg, 'avgoverchan', 'no');
cfg.latency     = ft_getopt(cfg, 'latency',     'all');
cfg.avgovertime = ft_getopt(cfg, 'avgovertime', 'no');

if isempty(cfg.parameter)
  if isfield(varargin{1}, 'trial')
    cfg.parameter = 'trial';
  elseif isfield(varargin{1}, 'individual')
    cfg.parameter = 'individual';
  elseif isfield(varargin{1}, 'avg')
    cfg.parameter = 'avg';
  end
end

% ensure that the data in all inputs has the same channels, time-axis, etc.
tmpcfg = keepfields(cfg, {'latency', 'avgovertime', 'channel', 'avgoverchan', 'parameter', 'showcallinfo'});
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

if strcmp(cfg.correctm, 'cluster') && length(varargin{1}.label)>1
  % this is required for clustering with multiple channels
  ft_checkconfig(cfg, 'required', 'neighbours');
end

dimord = getdimord(varargin{1}, cfg.parameter);
dimtok = tokenize(dimord, '_');
dimsiz = getdimsiz(varargin{1}, cfg.parameter);
dimsiz(end+1:length(dimtok)) = 1; % there can be additional trailing singleton dimensions
rptdim = find( strcmp(dimtok, 'subj') |  strcmp(dimtok, 'rpt') |  strcmp(dimtok, 'rpttap'));
datdim = find(~strcmp(dimtok, 'subj') & ~strcmp(dimtok, 'rpt') & ~strcmp(dimtok, 'rpttap'));
datsiz = dimsiz(datdim);

% pass these fields to the low-level functions, they should be removed further down
cfg.dimord = sprintf('%s_', dimtok{datdim});
cfg.dimord = cfg.dimord(1:end-1); % remove trailing _
cfg.dim    = dimsiz(datdim);

if isempty(rptdim)
  % repetitions are across multiple inputs
  dat = nan(prod(dimsiz), length(varargin));
  for i=1:length(varargin)
    tmp = varargin{i}.(cfg.parameter);
    dat(:,i) = tmp(:);
  end
else
  % repetitions are within inputs
  dat = cell(size(varargin));
  for i=1:length(varargin)
    tmp = varargin{i}.(cfg.parameter);
    if rptdim~=1
      % move the repetitions to the first dimension
      tmp = permute(tmp, [rptdim datdim]);
    end
    dat{i} = reshape(tmp, size(tmp,1), []);
  end
  dat = cat(1, dat{:});   % repetitions along 1st dimension
  dat = dat';             % repetitions along 2nd dimension
end

if size(cfg.design,2)~=size(dat,2)
  cfg.design = transpose(cfg.design);
end

design = cfg.design;

% determine the function handle to the intermediate-level statistics function
if exist(['ft_statistics_' cfg.method], 'file')
  statmethod = str2func(['ft_statistics_' cfg.method]);
else
  ft_error('could not find the corresponding function for cfg.method="%s"\n', cfg.method);
end
fprintf('using "%s" for the statistical testing\n', func2str(statmethod));

% check that the design completely describes the data
if size(dat,2) ~= size(cfg.design,2)
  ft_error('the length of the design matrix (%d) does not match the number of observations in the data (%d)', size(cfg.design,2), size(dat,2));
end

% determine the number of output arguments
try
  % the nargout function in MATLAB 6.5 and older does not work on function handles
  num = nargout(statmethod);
catch
  num = 1;
end

% perform the statistical test
if strcmp(func2str(statmethod),'ft_statistics_montecarlo')
  % because ft_statistics_montecarlo (or to be precise, clusterstat) requires to know whether it is getting source data,
  % the following (ugly) work around is necessary
  if num>1
    [stat, cfg] = statmethod(cfg, dat, design);
    cfg         = rollback_provenance(cfg); % ensure that changes to the cfg are passed back to the right level
  else
    [stat] = statmethod(cfg, dat, design);
  end
else
  if num>1
    [stat, cfg] = statmethod(cfg, dat, design);
    cfg         = rollback_provenance(cfg); % ensure that changes to the cfg are passed back to the right level
  else
    [stat] = statmethod(cfg, dat, design);
  end
end

if ~isstruct(stat)
  % only the probability was returned as a single matrix, reformat into a structure
  stat = struct('prob', stat);
end

% the statistical output contains multiple elements, e.g. F-value, beta-weights and probability
fn = fieldnames(stat);

for i=1:length(fn)
  if numel(stat.(fn{i}))==prod(datsiz)
    % reformat into the same dimensions as the input data
    stat.(fn{i}) = reshape(stat.(fn{i}), [datsiz 1]);
  end
end

% describe the dimensions of the output data
stat.dimord = cfg.dimord;

% copy the descripive fields into the output
stat = copyfields(varargin{1}, stat, {'time', 'label'});

% these were only present to inform the low-level functions
cfg = removefields(cfg, {'dim', 'dimord'});

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   varargin
ft_postamble provenance stat
ft_postamble history    stat
ft_postamble savevar    stat
