function stat = ft_statfun_bayesfactor(cfg, dat, design)

% FT_STATFUN_BAYESFACTOR computes the Bayes factor for a H0 of the data in two
% conditions having the same mean, versus H1 of the data having different means. This
% function supports both unpaired and paired designs and assumes flat priors.
%
% Lee and Wagenmakers (2013) provide these guidelines for its interpretation
%   IF B10 IS...    THEN YOU HAVE...
%     > 100           Extreme evidence for H1
%     30 – 100        Very strong evidence for H1
%     10 – 30         Strong evidence for H1
%     3 – 10          Moderate evidence for H1
%     1 – 3           Anecdotal evidence for H1
%     1               No evidence
%     1/3 – 1         Anecdotal evidence for H0
%     1/3 – 1/10      Moderate evidence for H0
%     1/10 – 1/30     Strong evidence for H0
%     1/30 – 1/100    Very strong evidence for H0
%     < 1/100         Extreme evidence for H0
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'ft_statfun_bayesfactor'
%
% The experimental design is specified as:
%   cfg.ivar  = independent variable, row number of the design that contains the labels of the conditions to be compared (default=1)
%   cfg.uvar  = optional, row number of design that contains the labels of the units-of-observation, i.e. subjects or trials (default=2)
%
% The labels for the independent variable should be specified as the number 1 and 2.
% The labels for the unit of observation should be integers ranging from 1 to the
% total number of observations (subjects or trials).
%
% The cfg.uvar option is only needed for paired data, you should leave it empty
% for non-paired data.
%
% See https://www.statisticshowto.datasciencecentral.com/bayes-factor-definition/ for some background.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS

% Copyright (C) 2020, Robert Oostenveld
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

% ensure that the required toolbox dependency is on the path
ft_hastoolbox('bayesFactor', 1);

% set the defaults
cfg.ivar = ft_getopt(cfg, 'ivar', 1);
cfg.uvar = ft_getopt(cfg, 'uvar', []); % default is empty, which means non-paired

if isempty(cfg.uvar)
  sel1 = find(design(cfg.ivar,:)==1); % select replications that belong to condition 1
  sel2 = find(design(cfg.ivar,:)==2); % select replications that belong to condition 2
  n1 = length(sel1);
  n2 = length(sel2);
  
  x1 = dat(:,sel1);
  x2 = dat(:,sel2);
  
  nchan = size(x1, 1);
  bf10  = nan(nchan, 1);
  prob  = nan(nchan, 1);
  ci_lo = nan(nchan, 1);
  ci_hi = nan(nchan, 1);
  tstat = nan(nchan, 1);
  diff  = mean(x1, 2) - mean(x2, 2); % this is the raw effect size
  
  for i=1:nchan
    [bf10(i), prob(i), CI, stats] = bf.ttest2(x1(i,:), x2(i,:));
    ci_lo(i) = CI(1);
    ci_hi(i) = CI(2);
    tstat(i) = stats.tstat;
  end
  
else
  subj = unique(design(cfg.uvar,:)); % it can also be paired over trials
  
  n = length(subj);
  sel1 = nan(size(subj));
  sel2 = nan(size(subj));
  for i=1:n
    sel1(i) = find(design(cfg.uvar,:)==subj(i) & design(cfg.ivar,:)==1);
    sel2(i) = find(design(cfg.uvar,:)==subj(i) & design(cfg.ivar,:)==2);
  end
  x1 = dat(:,sel1);
  x2 = dat(:,sel2);
  
  nchan = size(x1, 1);
  bf10  = nan(nchan, 1);
  prob  = nan(nchan, 1);
  ci_lo = nan(nchan, 1);
  ci_hi = nan(nchan, 1);
  tstat = nan(nchan, 1);
  diff  = mean(x1 - x2, 2); % this is the raw effect size
  
  for i=1:nchan
    [bf10(i), prob(i), CI, stats] = bf.ttest(x1(i,:) - x2(i,:));
    ci_lo(i) = CI(1);
    ci_hi(i) = CI(2);
    tstat(i) = stats.tstat;
  end
  
end % if uvar

% return the results for all channel-time-frequency points, for all voxels, or for all source locations
stat.bf10  = bf10(:);
stat.prob  = prob(:);
stat.diff  = diff(:);
stat.ci_lo = ci_lo(:);
stat.ci_hi = ci_hi(:);
stat.tstat = tstat(:);
