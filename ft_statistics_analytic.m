function [stat, cfg] = ft_statistics_analytic(cfg, dat, design)

% FT_STATISTICS_ANALYTIC performs a parametric statistical test on the
% data, based on a known (i.e. analytic) distribution of the test
% statistic. This function should not be called directly, instead
% you should call the function that is associated with the type of
% data on which you want to perform the test.
%
% Use as
%   stat = ft_timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = ft_freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = ft_sourcestatistics  (cfg, data1, data2, data3, ...)
% where the data is obtained from FT_TIMELOCKANALYSIS, FT_FREQANALYSIS
% or FT_SOURCEANALYSIS respectively, or from FT_TIMELOCKGRANDAVERAGE,
% FT_FREQGRANDAVERAGE or FT_SOURCEGRANDAVERAGE respectively.
%
% The configuration can contain
%   cfg.statistic        = string, statistic to compute for each sample or voxel (see below)
%   cfg.correctm         = apply multiple-comparison correction, 'no', 'bonferoni', 'holms', 'fdr' (default = 'no')
%   cfg.alpha            = critical value for rejecting the null-hypothesis (default = 0.05)
%   cfg.tail             = -1, 1 or 0 (default = 0)
%   cfg.ivar             = number or list with indices, independent variable(s)
%   cfg.uvar             = number or list with indices, unit variable(s)
%   cfg.wvar             = number or list with indices, within-block variable(s)
%
% The parametric statistic that is computed for each sample (and for
% which the analytic probability of the null-hypothesis is computed) is
% specified as
%   cfg.statistic       = 'indepsamplesT'     independent samples T-statistic,
%                         'indepsamplesF'     independent samples F-statistic,
%                         'indepsamplesregrT' independent samples regression coefficient T-statistic,
%                         'indepsamplesZcoh'  independent samples Z-statistic for coherence,
%                         'depsamplesT'       dependent samples T-statistic,
%                         'depsamplesF'       dependent samples F-statistic,
%                         'depsamplesregrT'   dependent samples regression coefficient T-statistic,
%                         'actvsblT'          activation versus baseline T-statistic.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS, FT_SOURCESTATISTICS

% Copyright (C) 2006, Robert Oostenveld
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

% set the defaults
if ~isfield(cfg, 'correctm'), cfg.correctm = 'no'; end
if ~isfield(cfg, 'alpha'),    cfg.alpha = 0.05;    end
if ~isfield(cfg, 'tail'),     cfg.tail = 0;        end

% determine the function handle to the low-level statistics function
if exist(['statistics_' cfg.statistic])
  statfun = str2func(['statistics_' cfg.statistic]);
elseif exist(['statfun_' cfg.statistic])
  statfun = str2func(['statfun_' cfg.statistic]);
else
  error(sprintf('could not find the statistics function "%s"\n', ['statfun_' cfg.statistic]));
end
fprintf('using ''%s'' to compute the single-sample statistic and probability\n', func2str(statfun));

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

switch lower(cfg.correctm)
  case 'bonferoni'
    fprintf('performing Bonferoni correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    stat.mask = stat.prob<=(cfg.alpha ./ numel(stat.prob));
  case 'holms'
    % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
    fprintf('performing Holms correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    [pvals, indx] = sort(stat.prob(:));                     % this sorts the significance probabilities from smallest to largest
    mask = pvals<=(cfg.alpha ./ ((length(pvals):-1:1)'));    % compare each significance probability against its individual threshold
    stat.mask = zeros(size(stat.prob));
    stat.mask(indx) = mask;
  case 'fdr'
    fprintf('performing FDR correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    stat.mask = fdr(stat.prob, cfg.alpha);
  otherwise
    fprintf('not performing a correction for multiple comparisons\n');
    stat.mask = stat.prob<=cfg.alpha;
end

