function [stat, cfg] = ft_statistics_analytic(cfg, dat, design)

% FT_STATISTICS_ANALYTIC performs a parametric statistical test on the data, based on
% a known (i.e. analytic) distribution of the test statistic. It is not recommended 
% to call this function directly, but instead you should call the function that is 
% associated with the type of data on  which you want to perform the test. This is 
% because important data bookkeeping  operations on the data are performed in the 
% higher-level functions, which are assumed to have been handled correctly for the 
% input arguments into this function. 
%
% Use as
%   stat = ft_timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = ft_freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = ft_sourcestatistics  (cfg, data1, data2, data3, ...)
%
% where the data is obtained from FT_TIMELOCKANALYSIS, FT_FREQANALYSIS or
% FT_SOURCEANALYSIS respectively, or from FT_TIMELOCKGRANDAVERAGE,
% FT_FREQGRANDAVERAGE or FT_SOURCEGRANDAVERAGE respectively 
% and with cfg.method = 'analytic'
%
% The configuration options that can be specified are:
%   cfg.statistic        = string, statistic to compute for each sample or voxel (see below)
%   cfg.correctm         = string, apply multiple-comparison correction, 'no', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
%   cfg.alpha            = number, critical value for rejecting the null-hypothesis (default = 0.05)
%   cfg.tail             = number, -1, 1 or 0 (default = 0)
%   cfg.ivar             = number or list with indices, independent variable(s)
%   cfg.uvar             = number or list with indices, unit variable(s)
%   cfg.wvar             = number or list with indices, within-block variable(s)
%
% The parametric statistic that is computed for each sample (and for
% which the analytic probability of the null-hypothesis is computed) is
% specified as
%   cfg.statistic       = 'indepsamplesT'           independent samples T-statistic,
%                         'indepsamplesF'           independent samples F-statistic,
%                         'indepsamplesregrT'       independent samples regression coefficient T-statistic,
%                         'indepsamplesZcoh'        independent samples Z-statistic for coherence,
%                         'depsamplesT'             dependent samples T-statistic,
%                         'depsamplesFmultivariate' dependent samples F-statistic MANOVA,
%                         'depsamplesregrT'         dependent samples regression coefficient T-statistic,
%                         'actvsblT'                activation versus baseline T-statistic.
% or you can specify your own low-level statistical function.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS, FT_SOURCESTATISTICS
% FT_STATISTICS_MONTECARLO, FT_STATISTICS_STATS, FT_STATISTICS_MVPA,
% FT_STATISTICS_CROSSVALIDATE

% Copyright (C) 2006-2015, Robert Oostenveld
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

% do a sanity check on the input data
assert(isnumeric(dat),    'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');
assert(isnumeric(design), 'this function requires numeric data as input, you probably want to use FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS instead');

% check whether the function has been called from ft_timelockstatistics, ft_freqstatistics, or ft_sourcestatistics
st = dbstack;
m  = mfilename;
if isscalar(st)
  ft_warning('It seems that %s has been called directly from the command line. This is not recommended, unless you know what you are doing', m);
elseif numel(st)>1 && ~ismember(st(2).name, {'ft_freqstatistics' 'ft_timelockstatistics' 'ft_sourcestatistics'})
  ft_warning('It seems that %s has not been called from one of the FT_XXXSTATISTICS functions. This is not recommended, unless you know what you are doing', m);
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval',  {'correctm', 'bonferoni', 'bonferroni'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'correctm', 'holms', 'holm'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'correctm', 'none', 'no'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'statfun', 'depsamplesF', 'ft_statfun_depsamplesFmultivariate'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'statfun', 'ft_statfun_depsamplesF', 'ft_statfun_depsamplesFmultivariate'});

% set the defaults
cfg.correctm = ft_getopt(cfg, 'correctm', 'no');
cfg.alpha    = ft_getopt(cfg, 'alpha', 0.05);
cfg.tail     = ft_getopt(cfg, 'tail', 0);

% fetch function handle to the low-level statistics function
statfun = ft_getuserfun(cfg.statistic, 'statfun');
if isempty(statfun)
  ft_error('could not locate the appropriate statistics function');
else
  ft_info('using "%s" for the single-sample statistics\n', func2str(statfun));
end

% tell the low-level statfun to compute the statistic, critical values, and the probability
cfg.computestat    = 'yes';
cfg.computecritval = 'yes';
cfg.computeprob    = 'yes';
% perform the statistical test
if nargout(statfun)>1
  [stat, cfg] = statfun(cfg, dat, design);
else
  [stat] = statfun(cfg, dat, design);
end
% these are only intenced for the low-level statfun, and should not be returned
cfg = rmfield(cfg, 'computestat');
cfg = rmfield(cfg, 'computeprob');
cfg = rmfield(cfg, 'computecritval');

if ~isfield(stat, 'prob')
  ft_warning('probability was not computed');
else
  switch lower(cfg.correctm)
    case 'bonferroni'
      ft_notice('performing Bonferoni correction for multiple comparisons\n');
      ft_notice('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
      stat.mask = stat.prob<=(cfg.alpha ./ numel(stat.prob));
    case 'holm'
      % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
      ft_notice('performing Holm-Bonferroni correction for multiple comparisons\n');
      ft_notice('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
      [pvals, indx] = sort(stat.prob(:));                                   % this sorts the significance probabilities from smallest to largest
      k = find(pvals > (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'first'); % compare each significance probability against its individual threshold
      mask = (1:length(pvals))'<k;
      stat.mask = zeros(size(stat.prob));
      stat.mask(indx) = mask;
    case 'hochberg'
      % test the most significant significance probability against alpha/N, the second largest against alpha/(N-1), etc.
      ft_notice('performing Hochberg''s correction for multiple comparisons (this is *not* the Benjamini-Hochberg FDR procedure!)\n');
      ft_notice('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
      [pvals, indx] = sort(stat.prob(:));                                   % this sorts the significance probabilities from smallest to largest
      k = find(pvals <= (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'last'); % compare each significance probability against its individual threshold
      mask = (1:length(pvals))'<=k;
      stat.mask = zeros(size(stat.prob));
      stat.mask(indx) = mask;
    case 'fdr'
      ft_notice('performing FDR correction for multiple comparisons\n');
      ft_notice('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
      stat.mask = fdr(stat.prob, cfg.alpha);
    case 'no'
      ft_notice('not performing a correction for multiple comparisons\n');
      stat.mask = stat.prob<=cfg.alpha;
    otherwise
      ft_error('unsupported option "%s" for cfg.correctm', cfg.correctm);
  end
end
