function [s,cfg] = statfun_indepsamplesT(cfg, dat, design);

% STATFUN_indepsamplesT calculates the independent samples T-statistic 
% on the biological data in dat (the dependent variable), using the information on 
% the independent variable (iv) in design.
%
% The external interface of this function has to be
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
% $Log: statfun_indepsamplesT.m,v $
% Revision 1.1  2008/09/23 07:31:15  roboos
% Moved all statfuns and trialfuns to their own directories, where they will be easier to find for the end-user. Also updated fieldtripdefs accordingly.
%
% Revision 1.10  2007/01/10 08:55:01  jansch
% fixed bug on line 97 and added a comment on why nreplc should be vectors
%
% Revision 1.9  2006/12/22 11:36:49  jansch
% *** empty log message ***
%
% Revision 1.8  2006/09/12 12:13:07  roboos
% removed default values for cfg.ivar and uvar, defaults should be specified elsewhere
%
% Revision 1.7  2006/06/07 12:51:18  roboos
% renamed cfg.ivrownr into cfg.ivar
% renamed cfg.uorownr into cfg.uvar
% renamed pval into prob for consistency with other fieldtrip functions
%
% Revision 1.6  2006/05/16 10:26:38  erimar
% Corrected bug involving sample sizes.
%
% Revision 1.5  2006/05/12 15:32:40  erimar
% Added functionality to calculate one- and two-sided critical values and
% p-values.
%
% Revision 1.4  2006/05/05 13:04:14  erimar
% Renamed several options and variables to make them consistent with other statfuns.
%
% Revision 1.3  2006/04/12 12:13:44  roboos
% added degrees of freedom to output if critvals is requested
%
% Revision 1.2  2006/04/11 16:17:27  roboos
% renamed some cfg options, changed the handling of whether to compute stat and/or critical values, changed the output into a stucture containing stat and/or critval
%

% set defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end;
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end;
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end;
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end;
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end;

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
  % probabilities can only be calculated if the test statistics are calculated
  cfg.computestat = 'yes';
end;

% perform some checks on the design
sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);
nreplc1 = sum(~isnan(dat(:,sel1)), 2);
nreplc2 = sum(~isnan(dat(:,sel2)), 2);
nrepl   = nreplc1 + nreplc2;
if any(nrepl<size(design,2)),
  warning('Not all replications are used for the computation of the statistic.');
end;
%if nrepl<3
%    error('The data must contain at least three trials/subjects.');
%end;
df = nrepl - 2;

if strcmp(cfg.computestat, 'yes')
  % compute the statistic
  avg1 = nanmean(dat(:,sel1), 2);
  avg2 = nanmean(dat(:,sel2), 2);
  var1 = nanstd(dat(:,sel1), 0, 2).^2;
  var2 = nanstd(dat(:,sel2), 0, 2).^2;
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
  if cfg.tail==-1
    s.critval = tinv(cfg.alpha,df);
  elseif  cfg.tail==0
    s.critval = [tinv(cfg.alpha/2,df),tinv(1-cfg.alpha/2,df)];
  elseif cfg.tail==1
    s.critval = tinv(1-cfg.alpha,df);
  end;
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
  end;
end
