function [s,cfg] = statfun_depsamplesregrT(cfg, dat, design);

% STATFUN_depsamplesregrT calculates dependent samples regression T-statistic 
% on the biological data in dat (the dependent variable), using the information on 
% the independent variable (iv) in design.
%
% The external interface of this function has to be
%   [s,cfg] = statfun_depsamplesregrT(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variable (iv) and the unit-of-observation (UO) 
%          factor,  Nreplications x Nvar
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
%   cfg.ivar        = row number of the design that contains the independent variable.
%   cfg.uvar        = row number of design that contains the labels of the UOs (subjects or trials)
%                        (default=2). The labels are assumed to be integers ranging from 1 to 
%                        the number of UOs.
%

% Copyright (C) 2006, Eric Maris
%
% $Log: statfun_depsamplesregrT.m,v $
% Revision 1.1  2008/09/23 07:31:15  roboos
% Moved all statfuns and trialfuns to their own directories, where they will be easier to find for the end-user. Also updated fieldtripdefs accordingly.
%
% Revision 1.10  2007/08/24 10:40:28  erimar
% Correct a bug (omission of a unit selection variable), detected by
% Vladimir Litvak,
%
% Revision 1.9  2007/08/16 10:15:31  erimar
% Added functionality with respect to control variables (specified in
% cfg.design and addressed via cfg.cvar). The contribution of control
% variables to the dependent variable (the biological data) can be
% partialled out.
%
% Revision 1.8  2007/08/16 10:12:30  erimar
% *** empty log message ***
%
% Revision 1.7  2006/09/12 12:13:07  roboos
% removed default values for cfg.ivar and uvar, defaults should be specified elsewhere
%
% Revision 1.6  2006/06/07 12:51:18  roboos
% renamed cfg.ivrownr into cfg.ivar
% renamed cfg.uorownr into cfg.uvar
% renamed pval into prob for consistency with other fieldtrip functions
%
% Revision 1.5  2006/05/17 11:59:55  erimar
% Corrected bugs after extensive checking of the properties of this
% statfun.
%
% Revision 1.4  2006/05/12 15:32:40  erimar
% Added functionality to calculate one- and two-sided critical values and
% p-values.
%
% Revision 1.3  2006/05/05 13:08:54  erimar
% Renamed several options and variables to make this statfun consistent with other
% statfuns.
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
if strcmp(cfg.computeprob,'yes') & strcmp(cfg.computestat,'no')
    error('P-values can only be calculated if the test statistics are calculated.');
end;

if ~isempty(cfg.cvar)
  condlabels=unique(design(cfg.cvar,:));
  nblocks=length(condlabels);
else
  nblocks=1;
end;

nunits = max(design(cfg.uvar,:));
df = nunits - 1;
if nunits<2
    error('The data must contain at least two units-of-observation (usually subjects).')
end;

if strcmp(cfg.computestat,'yes')
% compute the statistic
  regrweights=zeros(size(dat,1),nunits);
  for indx=1:nunits
    unitselvec=find(design(cfg.uvar,:)==indx);
    indvar=design(cfg.ivar,unitselvec); 
    if isempty(cfg.cvar)
      designmat=[ones(1,length(indvar));indvar];
    else
      designmat=zeros((nblocks+1),length(indvar));
      for blockindx=1:nblocks
        blockselvec=find(design(cfg.cvar,unitselved)==condlabels(blockindx));
        designmat(blockindx,blockselvec)=1;
      end;
      designmat((nblocks+1),:)=indvar;
    end;
    coeff=(designmat*designmat')\(designmat*dat(:,unitselvec)');
    regrweights(:,indx)=coeff((nblocks+1),:)';
  end;
  avgw=mean(regrweights,2);
  varw=var(regrweights,0,2);
  s.stat=sqrt(nunits)*avgw./sqrt(varw);
end;

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

