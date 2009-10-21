function [s] = statfun_indepsamplesZcoh(cfg, dat, design);

% STATFUN_indepsamplesZcoh calculates the independent samples coherence Z-statistic 
% on the biological data in dat (the dependent variable), using the information on 
% the independent variable (iv) in design.
%
% The external interface of this function has to be
%   [s,cfg] = statfun_indepsamplesT(cfg, dat, design);
% where
%   dat    contains the biological data, Nsamples x Nreplications
%          dat must contain fourier representations. 
%   design contains the independent variable (iv), Nreplications x Nvar
%
% The samples-dimension of the dat-variable must be the result of a
% reshaping-operation applied to a data structure with dimord
% chan_freq_time. The configuration must contain channel labels in
% cfg.label. This information is used to determine the number of channels. 
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
%

% Copyright (C) 2006, Eric Maris
%
% $Log: statfun_indepsamplesZcoh.m,v $
% Revision 1.1  2008/09/23 07:31:15  roboos
% Moved all statfuns and trialfuns to their own directories, where they will be easier to find for the end-user. Also updated fieldtripdefs accordingly.
%
% Revision 1.6  2006/09/12 12:13:07  roboos
% removed default values for cfg.ivar and uvar, defaults should be specified elsewhere
%
% Revision 1.5  2006/06/12 08:29:36  erimar
% Removed the code that added cfg.channelcmb to the configuration.
%
% Revision 1.4  2006/06/07 12:51:18  roboos
% renamed cfg.ivrownr into cfg.ivar
% renamed cfg.uorownr into cfg.uvar
% renamed pval into prob for consistency with other fieldtrip functions
%
% Revision 1.3  2006/05/17 11:59:55  erimar
% Corrected bugs after extensive checking of the properties of this
% statfun.
%
% Revision 1.2  2006/05/12 15:32:40  erimar
% Added functionality to calculate one- and two-sided critical values and
% p-values.
%
% Revision 1.1  2006/05/05 13:06:47  erimar
% First version of statfun_indepsamplesZcoh.
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

% perform some checks on the design
selc1 = find(design(cfg.ivar,:)==1);
selc2 = find(design(cfg.ivar,:)==2);
nreplc1 = length(selc1);
nreplc2 = length(selc2);
nrepl = nreplc1 + nreplc2;
if nrepl<size(design,1)
  error('Invalid specification of the independent variable in the design array.');
end;
if nreplc1<2 | nreplc2<2
    error('Every condition must contain at least two trials/tapers.');
end;
dfc1 = nreplc1*2;
dfc2 = nreplc2*2;

if strcmp(cfg.computestat, 'yes')
  % compute the statistic
  nsamples = size(dat,1);
  nchan = length(cfg.label);
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
end
