function [s, cfg] = ft_statfun_depsamplesregrT(cfg, dat, design)

% FT_STATFUN_DEPSAMPLESREGRT calculates dependent samples regression T-statistic 
% on the biological data in dat (the dependent variable), using the information on 
% the independent variable (ivar) in design.
%
% Use this function by calling one of the high-level statistics functions as
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option
%   cfg.statistic = 'ft_statfun_depsamplesregrT'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = ft_statfun_depsamplesregrT(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variable (ivar) and the unit-of-observation (uvar) 
%          factor,  Nreplications x Nvar
%
% Configuration options
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
%
% The following options are relevant if cfg.computecritval='yes' and/or
% cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail  = -1, 0, or 1, left, two-sided, or right (default=1)
%               cfg.tail in combination with cfg.computecritval='yes'
%               determines whether the critical value is computed at
%               quantile cfg.alpha (with cfg.tail=-1), at quantiles
%               cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%               quantile (1-cfg.alpha) (with cfg.tail=1).
%
% Design specification
%   cfg.ivar  = row number of the design that contains the independent variable.
%   cfg.uvar  = row number of design that contains the labels of the units-of-observation (subjects or trials)
%               (default=2). The labels are assumed to be integers ranging from 1 to 
%               the number of units-of-observation.

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
if ~isfield(cfg,'uvar') || isempty(cfg.uvar)
    error('uvar must be specified for dependent samples statistics');
end

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

