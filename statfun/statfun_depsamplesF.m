function [s,cfg] = statfun_depsamplesF(cfg, dat, design)

% STATFUN_depsamplesF calculates the dependent samples F-statistic 
% on the biological data in dat (the dependent variable), using the information on 
% the independent variable (iv) in design.
%
% Use this function by calling one of the high-level statistics functions as:
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'depsamplesF'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = statfun_depsamplesF(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%   design contains the independent variable (iv) and the unit-of-observation (UO) 
%          factor,  Nfac x Nreplications
%
% Configuration options:
%   cfg.contrastcoefs  = matrix of contrast coefficients determining the
%                      effect being tested. The number of columns of this
%                      matrix has to be equal to the number of conditions. 
%                      The default is a matrix that specifies the
%                      main effect of the independent variable. This matrix
%                      has size [(ncond-1),ncond]. 
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
%                        compared (default=1). The labels range from 1 to the number of conditions.
%   cfg.uvar        = row number of design that contains the labels of the UOs (subjects or trials)
%                        (default=2). The labels are assumed to be integers ranging from 1 to 
%                        the number of UOs.
%

% Copyright (C) 2006, Eric Maris
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% set defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end;
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end;
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end;
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end;
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end;

nconds=length(unique(design(cfg.ivar,:)));
if ~isfield(cfg,'contrastcoefs')
    % specify the default contrast coefficient matrix.
    ncontrasts = nconds-1;
    cfg.contrastcoefs = zeros(ncontrasts,nconds);
    cfg.contrastcoefs(:,1) = 1;
    for contrastindx=1:ncontrasts
        cfg.contrastcoefs(contrastindx,contrastindx+1)=-1;
    end;
else
    ncontrasts = size(cfg.contrastcoefs,1);
end;

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') & strcmp(cfg.computestat,'no')
    error('P-values can only be calculated if the test statistics are calculated.');
end;
if ~isfield(cfg,'uvar') || isempty(cfg.uvar)
    error('uvar must be specified for dependent samples statistics');
end

% perform some checks on the design
nuospercond=zeros(nconds,1);
for condindx=1:nconds
    nuospercond(condindx)=length(find(design(cfg.ivar,:)==condindx));
end;
if sum(nuospercond)<size(design,2) | nuospercond~=(nuospercond(1)*ones(nconds,1))
  error('Invalid specification of the design array.');
end;
nunits = max(design(cfg.uvar,:));
dfdenom = nunits - ncontrasts;
if dfdenom<1
    error('The data must contain more units-of-observation (usually subjects) than the number of contrasts.')
end;
nrepl=nunits*nconds;
if (nrepl~=sum(nuospercond)) | (nrepl~=size(dat,2))
  error('Invalid specification of the design array.');
end;
nsmpls = size(dat,1);

if strcmp(cfg.computestat,'yes')
% compute the statistic
    % store the positions of the condition labels nunits-by-nconds array
    poslabelsperunit=zeros(nunits,nconds);
    for condindx=1:nconds
        poslabel=find(design(cfg.ivar,:)==condindx);
        [dum,i]=sort(design(cfg.uvar,poslabel),'ascend');
        poslabelsperunit(:,condindx)=poslabel(i);
    end;
    % reshape poslabelsperunit into a row vector that contains the
    % replications of the first condition on the first nunits positions,
    % the replications of the second condition on the second nunits
    % positions, etc.
    poslabelsperunit=reshape(poslabelsperunit,1,nrepl);
    s.stat=zeros(nsmpls,1);
    for smplindx=1:nsmpls
        datonesmpl=reshape(dat(smplindx,poslabelsperunit),nunits,nconds);
        contrasts=datonesmpl*cfg.contrastcoefs';
        contrastavg=mean(contrasts,1);
        dev=contrasts-repmat(contrastavg,nunits,1);
        covmat=(dev'*dev)/(nunits-1);
        s.stat(smplindx)=nunits*contrastavg*inv(covmat)*contrastavg';
    end;
end;

if strcmp(cfg.computecritval,'yes')
  % also compute the critical values
  s.dfnum = ncontrasts;
  s.dfdenom = nunits - ncontrasts;
  if cfg.tail==-1
      error('For a dependent samples F-statistic, it does not make sense to calculate a left tail critical value.');
  end;
  if cfg.tail==0
      error('For a dependent samples F-statistic, it does not make sense to calculate a two-sided critical value.');
  end;
  if cfg.tail==1
    s.critval = ((nunits-1).*ncontrasts./(nunits-ncontrasts)).*finv(1-cfg.alpha,s.dfnum,s.dfdenom);
  end;
end

if strcmp(cfg.computeprob,'yes')
  % also compute the p-values
  s.dfnum = ncontrasts;
  s.dfdenom = nunits - ncontrasts;
  if cfg.tail==-1
      error('For a dependent samples F-statistic, it does not make sense to calculate a left tail p-value.');
  end;
  if cfg.tail==0
      error('For a dependent samples F-statistic, it does not make sense to calculate a two-sided p-value.');
  end;
  if cfg.tail==1
    scaledstat = ((nunits-ncontrasts)./((nunits-1).*ncontrasts)).*s.stat;
    s.prob = 1-fcdf(scaledstat,s.dfnum,s.dfdenom);
  end;
end
