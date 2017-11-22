function [s, cfg] = ft_statfun_correlationT(cfg, dat, design)

% FT_STATFUN_CORRELATIONT calculates correlation coefficient T-statistics
% on the biological data in dat (the dependent variable), using the 
% information on the independent variable (predictor) in design. The
% correlation coefficients are stored in the rho field of output s.
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option
%   cfg.statistic = 'ft_statfun_correlationT'
%
% See FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = ft_statfun_correlationT(cfg, dat, design);
% where
%   dat    contains the biological data,  Nsamples x Nreplications
%   design contains the independent variable,  Nvar X Nreplications
%
% Configuration options
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
% The following options are relevant if cfg.computecritval='yes' and/or
% cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail  = -1, 0, or 1, left, two-sided, or right (default=1)
%               cfg.tail in combination with cfg.computecritval='yes'
%               determines whether the critical value is computed at
%               quantile cfg.alpha (with cfg.tail=-1), at quantiles
%               cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%               quantile (1-cfg.alpha) (with cfg.tail=1).
%   cfg.type  = 'Pearson' to compute Pearson's correlation (default), see 'help corr' for other options.
%
% Design specification
%   cfg.ivar  = row number of the design that contains the independent variable (default=1)

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

% set defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end
if ~isfield(cfg, 'type'),              cfg.type           = 'Pearson'; end

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
    ft_error('P-values can only be calculated if the test statistics are calculated.');
end
if isfield(cfg,'uvar') && ~isempty(cfg.uvar)
    ft_error('cfg.uvar should not exist for a correlation statistic');
end

[nsmpl,nrepl] = size(dat);
df = nrepl - 1;
if df<1
  ft_error('Insufficient error degrees of freedom for this analysis.')
end

if strcmp(cfg.computestat,'yes') % compute the statistic    
    % calculate the correlation coefficient between the dependent variable and the predictor
    rho = corr(dat', design', 'type', cfg.type);
    clear dat
    
    % convert correlation coefficient to t-statistic (for MCP correction): t^2 = DF*R^2 / (1-R^2)
    tstat = rho*(sqrt(nrepl-2))./sqrt((1-rho.^2));
    
    s.stat = tstat; % store t values in s.stat variable for use with ft_statistics_montecarlo.m
    s.rho = rho; % store r values in s.rho variable (these are the actual correlation coefficients)
    clear rho tstat
end

if strcmp(cfg.computecritval,'yes')
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

if strcmp(cfg.computeprob,'yes')
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
