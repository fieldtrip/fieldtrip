function [s, cfg] = ft_statfun_correlationT(cfg, dat, design)

% FT_STATFUN_CORRELATIONT calculates correlation coefficient T-statistics on the
% biological data in dat (the dependent variable), using the information on the
% independent variable (predictor) in design. The correlation coefficients themselves
% are stored in the output structure as rho.
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'ft_statfun_correlationT'
%
% You can specify the following configuration options:
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
%
% The following options are relevant if cfg.computecritval='yes' and/or cfg.computeprob='yes':
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail  = -1, 0, or 1, left, two-sided, or right (default=1)
%               cfg.tail in combination with cfg.computecritval='yes'
%               determines whether the critical value is computed at
%               quantile cfg.alpha (with cfg.tail=-1), at quantiles
%               cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%               quantile (1-cfg.alpha) (with cfg.tail=1)
%   cfg.type  = 'Pearson' to compute Pearson's correlation (default),
%                see 'help corr' for other options.
%
% The experimental design is specified as:
%   cfg.ivar  = row number of the design that contains the independent variable, i.e. the predictor (default=1)
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS

% Copyright (C) 2015, Arjen Stolk
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

% set the defaults
cfg.computestat    = ft_getopt(cfg, 'computestat', 'yes');
cfg.computecritval = ft_getopt(cfg, 'computecritval', 'no');
cfg.computeprob    = ft_getopt(cfg, 'computeprob', 'no');
cfg.alpha          = ft_getopt(cfg, 'alpha', 0.05);
cfg.tail           = ft_getopt(cfg, 'tail', 1);
cfg.type           = ft_getopt(cfg, 'type', 'Pearson'); % 'Pearson', 'Kendall' or 'Spearman'
cfg.ivar           = ft_getopt(cfg, 'ivar', 1);

% perform some checks on the configuration
if strcmp(cfg.computeprob, 'yes') && strcmp(cfg.computestat, 'no')
  ft_error('P-values can only be calculated if the test statistics are calculated.');
end

[nsmpl, nrepl] = size(dat);
df = nrepl - 1;
if df<1
  ft_error('Insufficient degrees of freedom for this analysis.')
end

if strcmp(cfg.computestat, 'yes') % compute the statistic
  % calculate the correlation coefficient between the dependent variable and the predictor
  rho = corr(dat', design', 'type', cfg.type);

  % convert correlation coefficient to t-statistic (for MCP correction): t^2 = DF*R^2 / (1-R^2)
  tstat = rho*(sqrt(nrepl-2))./sqrt((1-rho.^2));

  s.stat = tstat; % store t values in s.stat variable, these are used for the inference in FT_STATISTICS_MONTECARLO
  s.rho = rho;    % store r values in s.rho variable, these are the actual correlation coefficients
end

if strcmp(cfg.computecritval, 'yes')
  % also compute the critical values
  s.df      = df;
  if cfg.tail==-1
    s.critval = tinv(cfg.alpha,df);
  elseif  cfg.tail==0
    s.critval = [tinv(cfg.alpha/2,df),tinv(1-cfg.alpha/2,df)];
  elseif cfg.tail==1
    s.critval = tinv(1-cfg.alpha,df);
  end
end

if strcmp(cfg.computeprob, 'yes')
  % also compute the p-values
  s.df      = df;
  if cfg.tail==-1
    s.prob = tcdf(s.stat,s.df);
  elseif  cfg.tail==0
    s.prob = 2*tcdf(-abs(s.stat),s.df);
  elseif cfg.tail==1
    s.prob = 1-tcdf(s.stat,s.df);
  end
end
