function [s, cfg] = statfun_indepsamplesT(cfg, dat, design)

% STATFUN_indepsamplesT calculates the independent samples T-statistic 
% on the biological data in dat (the dependent variable), using the information on 
% the independent variable (iv) in design.
%
% Use this function by calling one of the high-level statistics functions as:
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'indepsamplesT'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = statfun_indepsamplesT(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variable (iv),  Nfac x Nreplications
%
% Configuration options:
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
%
% The following options are relevant if cfg.computecritval='yes' and/or
% cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail = -1, 0, or 1, left, two-sided, or right (default=1)
%              cfg.tail in combination with cfg.computecritval='yes'
%              determines whether the critical value is computed at
%              quantile cfg.alpha (with cfg.tail=-1), at quantiles
%              cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%              quantile (1-cfg.alpha) (with cfg.tail=1).
%
% Design specification:
%   cfg.ivar        = row number of the design that contains the labels of the conditions that must be 
%                        compared (default=1). The labels are the numbers 1 and 2.
%

% Copyright (C) 2006, Eric Maris
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set the defaults
if ~isfield(cfg, 'computestat'),    cfg.computestat    = 'yes'; end
if ~isfield(cfg, 'computecritval'), cfg.computecritval = 'no';  end
if ~isfield(cfg, 'computeprob'),    cfg.computeprob    = 'no';  end
if ~isfield(cfg, 'alpha'),          cfg.alpha          = 0.05;  end
if ~isfield(cfg, 'tail'),           cfg.tail           = 1;     end

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
  % probabilities can only be calculated if the test statistics are calculated
  cfg.computestat = 'yes';
end;
if isfield(cfg,'uvar') && ~isempty(cfg.uvar)
    error('cfg.uvar should not exist for an independent samples statistic');
end

% perform some checks on the design
sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);
nreplc1 = sum(~isnan(dat(:,sel1)), 2);
nreplc2 = sum(~isnan(dat(:,sel2)), 2);
nrepl   = nreplc1 + nreplc2;
if any(nrepl<size(design,2)),
  warning_once('Not all replications are used for the computation of the statistic.');
end;
%if nrepl<3
%    error('The data must contain at least three trials/subjects.');
%end;
df = nrepl - 2;

if strcmp(cfg.computestat, 'yes')
  % compute the statistic use nanmean only if necessary
  if sum(any(isnan(dat(:,sel1))))
    avg1 = nanmean(dat(:,sel1), 2);
    var1 = nanstd(dat(:,sel1), 0, 2).^2;
  else
    avg1 = mean(dat(:,sel1), 2);
    var1 = std(dat(:,sel1), 0, 2).^2;
  end
  if sum(any(isnan(dat(:,sel2))))
    avg2 = nanmean(dat(:,sel2), 2);
    var2 = nanstd(dat(:,sel2), 0, 2).^2;
  else
    avg2 = mean(dat(:,sel2), 2);
    var2 = std(dat(:,sel2), 0, 2).^2;
  end
 
  varc = (1./nreplc1 + 1./nreplc2).*((nreplc1-1).*var1 + (nreplc2-1).*var2)./df;
  
  %in the case of non-equal triallengths, and tfrs as input-data nreplc are
  %vectors with different values. when the triallengths are equal, and the
  %input is a tfr, nreplc are vectors with either zeros (all trials contain nan
  %meaning that t_ftimwin did not fit around data), or the number of trials
  s.stat = (avg1 - avg2)./sqrt(varc);
end

if strcmp(cfg.computecritval,'yes')
  % also compute the critical values
  s.df      = df;
  if all(df==df(1)), df = df(1); end
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
