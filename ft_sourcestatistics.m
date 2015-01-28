function [stat] = ft_sourcestatistics(cfg, varargin)

% FT_SOURCESTATISTICS computes the probability for a given null-hypothesis using
% a parametric statistical test or using a non-parametric randomization test.
%
% Use as
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% where the input data is the result from FT_SOURCEANALYSIS, FT_SOURCEDESCRIPTIVES
% or FT_SOURCEGRANDAVERAGE.  The source structures should be spatially alligned
% to each other and should have the same positions for the source grid.
%
% The configuration should contain the following option for data selection
%   cfg.parameter  = string, describing the functional data to be processed, e.g. 'pow', 'nai' or 'coh'
%
% Furthermore, the configuration should contain:
%   cfg.method       = different methods for calculating the probability of the null-hypothesis,
%                    'montecarlo'    uses a non-parametric randomization test to get a Monte-Carlo estimate of the probability,
%                    'analytic'      uses a parametric test that results in analytic probability,
%                    'stats'         (soon deprecated) uses a parametric test from the MATLAB statistics toolbox,
%
% The other cfg options depend on the method that you select. You
% should read the help of the respective subfunction FT_STATISTICS_XXX
% for the corresponding configuration options and for a detailed
% explanation of each method.
%
% See also FT_SOURCEANALYSIS, FT_SOURCEDESCRIPTIVES, FT_SOURCEGRANDAVERAGE

% Deprecated cfg.method options:
%                    'parametric'    uses the MATLAB statistics toolbox (very similar to 'stats'),
%                    'randomization' uses randomization of the data prior to source reconstruction,
%                    'randcluster'   uses randomization of the data prior to source reconstruction
%                                    in combination with spatial clusters.

% FIXME the following needs to be reimplemented
%
% You can restrict the statistical analysis to regions of interest (ROIs)
% or to the average value inside ROIs using the following options:
%   cfg.atlas        = filename of the atlas
%   cfg.roi          = string or cell of strings, region(s) of interest from anatomical atlas
%   cfg.avgoverroi   = 'yes' or 'no' (default = 'no')
%   cfg.hemisphere   = 'left', 'right', 'both', 'combined', specifying this is
%                      required when averaging over regions
%   cfg.inputcoord   = 'mni' or 'tal', the coordinate system in which your source
%                      reconstruction is expressed

% Copyright (C) 2005-2014, Robert Oostenveld
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
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar varargin

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required',    {'method', 'design'});
cfg = ft_checkconfig(cfg, 'renamed',     {'approach',   'method'});


%%%%%%%%%

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required',    {'method', 'design'});
cfg = ft_checkconfig(cfg, 'renamed',     {'approach',   'method'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'transform'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'trials'}); % this used to be present until 24 Dec 2014, but was deemed too confusing by Robert

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'source', 'feedback', 'no');
  varargin{i} = ft_datatype_source(varargin{i}, 'version', 'upcoming');
end

% set the defaults
cfg.parameter   = ft_getopt(cfg, 'parameter'); % default is set below
cfg.correctm    = ft_getopt(cfg, 'correctm');
cfg.latency     = ft_getopt(cfg, 'latency',     'all');
cfg.avgovertime = ft_getopt(cfg, 'avgovertime', 'no');
cfg.frequency   = ft_getopt(cfg, 'frequency',   'all');
cfg.avgoverfreq = ft_getopt(cfg, 'avgoverfreq', 'no');

if isempty(cfg.parameter)
  if isfield(varargin{1}, 'pow')
    cfg.parameter = 'pow';
  end
end

if length(cfg.parameter)>4 && strcmp(cfg.parameter(1:4), 'avg.')
  cfg.parameter = cfg.parameter(5:end); % remove the 'avg.' part
end

% ensure that the data in all inputs has the same channels, time-axis, etc.
tmpcfg = keepfields(cfg, {'frequency', 'avgoverfreq', 'latency', 'avgovertime', 'parameter'});
[varargin{:}] = ft_selectdata(tmpcfg, varargin{:});
% restore the provenance information
[cfg, varargin{:}] = rollback_provenance(cfg, varargin{:});

dimord = getdimord(varargin{1}, cfg.parameter);
dimsiz = getdimsiz(varargin{1}, cfg.parameter);
dimtok = tokenize(dimord, '_');
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
if exist(['ft_statistics_' cfg.method])
  statmethod = str2func(['ft_statistics_' cfg.method]);
else
  error('could not find the corresponding function for cfg.method="%s"\n', cfg.method);
end
fprintf('using "%s" for the statistical testing\n', func2str(statmethod));

% check that the design completely describes the data
if size(dat,2) ~= size(cfg.design,2)
  error('the length of the design matrix (%d) does not match the number of observations in the data (%d)', size(cfg.design,2), size(dat,2));
end

% determine the number of output arguments
try
  % the nargout function in Matlab 6.5 and older does not work on function handles
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
  else
    [stat] = statmethod(cfg, dat, design);
  end
else
  if num>1
    [stat, cfg] = statmethod(cfg, dat, design);
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
stat = copyfields(varargin{1}, stat, {'freq', 'time', 'pos', 'dim', 'transform'});

% these were only present to inform the low-level functions
cfg = rmfield(cfg, 'dim');
cfg = rmfield(cfg, 'dimord');

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous varargin
ft_postamble history stat
ft_postamble savevar stat

