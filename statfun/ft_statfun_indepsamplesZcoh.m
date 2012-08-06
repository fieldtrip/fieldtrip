function [s, cfg] = statfun_indepsamplesZcoh(cfg, dat, design)

% STATFUN_indepsamplesZcoh calculates the independent samples coherence Z-statistic 
% on the biological data in dat (the dependent variable), using the information on 
% the independent variable (iv) in design.
%
% Use this function by calling one of the high-level statistics functions as:
%   [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%   [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
% with the following configuration option:
%   cfg.statistic = 'indepsamplesZcoh'
% see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = statfun_indepsamplesZcoh(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%          dat must contain fourier representations. 
%   design contains the independent variable (iv), Nreplications x Nvar
%
% The samples-dimension of the dat-variable must be the result of a
% reshaping-operation applied to a data structure with dimord
% chan_(freq_time) or pos_(freq_time). The configuration must contain
% channel labels in cfg.label or position information in cfg.pos. This
% information is used to determine the number of channels. 
% The dimord of the output fields is [prod(nchancmb,nfreq,ntime),1]. The
% channel combinations are the elements of the lower diagonal of the
% cross-spectral density matrix.
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
%   cfg.ivar        = column number of the design that contains the labels of the conditions that must be 
%                     compared (default=1). The labels are the numbers 1 and 2.

% Copyright (C) 2006, Eric Maris

% set the defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end;
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end;
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end;
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end;
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end;

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
    error('P-values can only be calculated if the test statistics are calculated.');
end;
if isfield(cfg,'uvar') && ~isempty(cfg.uvar)
    error('cfg.uvar should not exist for an independent samples statistic');
end
% if ~isfield(cfg, 'label') && ~isfield(cfg, 'pos')
%   error('the configuration needs to contain either a label or a pos field');
% elseif isfield(cfg, 'label') && isfield(cfg, 'pos') && ~isempty(cfg.label) && ~isempty(cfg.pos)
%   error('the configuration needs to contain either a non-empty label or a non-empty pos field');
% elseif isfield(cfg, 'label') && ~isempty(cfg.label)
%   nchan = length(cfg.label);
% elseif isfield(cfg, 'pos') && ~isempty(cfg.pos)
%   nchan = size(cfg.pos,1);
% end
nchan = cfg.dim(1);

% perform some checks on the design
selc1 = find(design(cfg.ivar,:)==1);
selc2 = find(design(cfg.ivar,:)==2);
nreplc1 = length(selc1);
nreplc2 = length(selc2);
nrepl = nreplc1 + nreplc2;
if nrepl<size(design,1)
  error('Invalid specification of the independent variable in the design array.');
end;
if nreplc1<2 || nreplc2<2
    error('Every condition must contain at least two trials/tapers.');
end;
dfc1 = nreplc1*2;
dfc2 = nreplc2*2;

if strcmp(cfg.computestat, 'yes')
  % compute the statistic
  nsamples = size(dat,1);
  %nchan = length(cfg.label); %this is computed earlier
  chancmbsel = find(tril(ones(nchan),-1));
  nfreqtim = nsamples/nchan;
  nchancmb = length(chancmbsel);
  nnewsamples = nchancmb*nfreqtim;
  s.stat = zeros(nnewsamples,1);
  for freqtimindx=1:nfreqtim
    chansel=((freqtimindx-1)*nchan + 1):(freqtimindx*nchan);
    csdc1=dat(chansel,selc1)*dat(chansel,selc1)'/nreplc1;
    powerc1=diag(csdc1);
    normmat=diag(powerc1.^(-1/2));
    csdc1=normmat*csdc1*normmat;
    csdc2=dat(chansel,selc2)*dat(chansel,selc2)'/nreplc2;
    powerc2=diag(csdc2);
    normmat=diag(powerc2.^(-1/2));
    csdc2=normmat*csdc2*normmat;
    biasc1=1./(dfc1-2);
    biasc2=1./(dfc2-2);
    denomZ=sqrt(1./(dfc1-2) + 1./(dfc2-2));
    tempstat=(atanh(abs(csdc1))-biasc1-atanh(abs(csdc2))+biasc2)./denomZ;
    s.stat(((freqtimindx-1)*nchancmb + 1):(freqtimindx*nchancmb))=tempstat(chancmbsel);
  end;
  s.stat = reshape(s.stat, [nchancmb nfreqtim]);
end;

if strcmp(cfg.computecritval,'yes')
  % also compute the critical values
  if cfg.tail==-1
    s.critval = norminv(cfg.alpha);
  elseif  cfg.tail==0
    s.critval = [norminv(cfg.alpha/2),norminv(1-cfg.alpha/2)];
  elseif cfg.tail==1
    s.critval = norminv(1-cfg.alpha);
  end;
end

if strcmp(cfg.computeprob,'yes')
  % also compute the p-values
  if cfg.tail==-1
    s.prob = normcdf(s.stat);
  elseif  cfg.tail==0
    s.prob = 2*normcdf(-abs(s.stat));
  elseif cfg.tail==1
    s.prob = 1-normcdf(s.stat);
  end;
  s.prob = reshape(s.prob, [nchancmb nfreqtim]);
end

% adjust the dimord
if strcmp(cfg.dimord(1:3), 'pos')
  cfg.dimord = ['poscmb',cfg.dimord(4:end)];
elseif strcmp(cfg.dimord(1:4), 'chan')
  cfg.dimord = ['chancmb',cfg.dimord(5:end)];
end

% append an indexing matrix to the cfg to be able to recover the channel
% combinations
chanindx = tril(true(nchan),-1);
cmbindx1 = repmat((1:nchan)', [1 nchan]);
cmbindx2 = repmat((1:nchan),  [nchan 1]);
cfg.chancmbindx(:,1) = cmbindx1(chanindx);
cfg.chancmbindx(:,2) = cmbindx2(chanindx);
