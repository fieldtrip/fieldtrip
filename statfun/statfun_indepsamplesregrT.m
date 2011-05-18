function [s,cfg] = statfun_indepsamplesregrT(cfg, dat, design)

% STATFUN_indepsamplesregrT calculates independent samples regression
% coefficient T-statistics on the biological data in dat (the dependent
% variable), using the information on the independent variable
% (predictor) in design.
%
% Use this function by calling one of the high-level statistics functions as:
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'indepsamplesregrT'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = statfun_indepsamplesregrT(cfg, dat, design);
% where
%   dat    contains the biological data,  Nsamples x Nreplications
%   design contains the independent variable,  Nreplications x Nvar
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
%   cfg.ivar        = row number of the design that contains the independent variable (default=1)

% Copyright (C) 2006, Eric Maris

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
if isfield(cfg,'uvar') && ~isempty(cfg.uvar)
    error('cfg.uvar should not exist for an independent samples statistic');
end

if ~isempty(cfg.cvar)
  condlabels=unique(design(cfg.cvar,:));
  nblocks=length(condlabels);
else
  nblocks=1;
end;

[nsmpl,nrepl] = size(dat);
df = nrepl - nblocks - 1;
if df<1
  error('Insufficient error degrees of freedom for this analysis.')
end;

if strcmp(cfg.computestat, 'yes')
  % compute the statistic
    indvar = design(cfg.ivar,:);
    if isempty(cfg.cvar)
      designmat = [ones(nrepl,1) indvar']; % designmat is a matrix of order Nrepl x 2
    else
      designmat = zeros(nrepl,(nblocks+1));
      for blockindx=1:nblocks
        selvec=find(design(cfg.cvar,:)==condlabels(blockindx));
        designmat(selvec,blockindx)=1; 
      end;
      designmat(:,(nblocks+1))=indvar'; % designmat is a matrix of order Nrepl x (nblocks+1)
    end;
    cpmat = designmat'*designmat;
    invcpmat = inv(cpmat);
    projmat = invcpmat*designmat';
    B = dat*projmat'; % B is a matrix of order Nsamples x (nblocks+1)
    res = dat - B*designmat';
    resvar = zeros(nsmpl,1);
    for indx=1:nsmpl
      resvar(indx)=res(indx,:)*res(indx,:)';
    end;
    resvar=resvar/df;
    
    se=sqrt(invcpmat(nblocks+1,nblocks+1)*resvar);
    s.stat=B(:,nblocks+1)./se;
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

