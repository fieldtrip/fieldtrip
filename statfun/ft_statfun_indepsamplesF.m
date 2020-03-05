function [s, cfg] = ft_statfun_indepsamplesF(cfg, dat, design)

% FT_STATFUN_INDEPSAMPLESF calculates the independent samples F-statistic on the
% biological data in dat (the dependent variable), using the information on the
% independent variable (ivar) in design.
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'ft_statfun_indepsamplesF'
%
% You can specify the following configuration options:
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default= 'yes')
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
%
% The experimental design is specified as:
%   cfg.ivar  = independent variable, row number of the design that contains the labels of the conditions to be compared (default=1)
%
% The labels for the independent variable should be specified as numbers ranging
% from 1 to the number of conditions.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS

% Copyright (C) 2006, Eric Maris
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
cfg.ivar           = ft_getopt(cfg, 'ivar', 1);

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
  ft_error('P-values can only be calculated if the test statistics are calculated.');
end
if isfield(cfg,'uvar') && ~isempty(cfg.uvar)
  ft_error('cfg.uvar should not exist for an independent samples statistic');
end

ncond = length(unique(design(cfg.ivar,:)));
nrepl = 0;
for condindx=1:ncond
  nrepl=nrepl+length(find(design(cfg.ivar,:)==condindx));
end
if nrepl<size(design,2)
  ft_error('Invalid specification of the independent variable in the design array.');
end
if nrepl<=ncond
  ft_error('The must be more trials/subjects than levels of the independent variable.');
end
dfnum   = ncond - 1;
dfdenom = nrepl - ncond;

nsmpls = size(dat,1);

if strcmp(cfg.computestat, 'yes')
  % compute the statistic
  nobspercell = zeros(1,ncond);
  avgs = zeros(nsmpls,ncond);
  pooledvar = zeros(nsmpls,1);
  for condindx=1:ncond
    sel=find(design(cfg.ivar,:)==condindx);
    nobspercell(condindx)=length(sel);
    avgs(:,condindx)=mean(dat(:,sel),2);
    pooledvar = pooledvar + nobspercell(condindx)*var(dat(:,sel),1,2);
  end
  pooledvar = pooledvar/dfdenom;
  globalavg = mean(dat,2);
  mseffect = ((avgs-repmat(globalavg,1,ncond)).^2)*nobspercell'./dfnum;
  s.stat = mseffect./pooledvar;
end

if strcmp(cfg.computecritval,'yes')
  % also compute the critical values
  s.dfnum   = dfnum;
  s.dfdenom = dfdenom;
  if cfg.tail==-1
    ft_error('For an independent samples F-statistic, it does not make sense to calculate a left tail critical value.');
  end
  if cfg.tail==0
    ft_error('For an independent samples F-statistic, it does not make sense to calculate a two-sided critical value.');
  end
  if cfg.tail==1
    s.critval = finv(1-cfg.alpha,s.dfnum,s.dfdenom);
  end
end

if strcmp(cfg.computeprob,'yes')
  % also compute the p-values
  s.dfnum   = dfnum;
  s.dfdenom = dfdenom;
  if cfg.tail==-1
    ft_error('For an independent samples F-statistic, it does not make sense to calculate a left tail p-value.');
  end
  if cfg.tail==0
    ft_error('For an independent samples F-statistic, it does not make sense to calculate a two-sided p-value.');
  end
  if cfg.tail==1
    s.prob = 1-fcdf(s.stat,s.dfnum,s.dfdenom);
  end
end
