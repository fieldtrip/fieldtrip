function p = prior_sqrtt(varargin)
%PRIOR_SQRTT  Student-t prior structure for the square root of the parameter
%       
%  Description
%    P = PRIOR_SQRTT('PARAM1', VALUE1, 'PARAM2', VALUE2, ...) 
%    creates for the square root of the parameter Student's
%    t-distribution prior structure in which the named parameters
%    have the specified values. Any unspecified parameters are set
%    to default values.
%
%    P = PRIOR_SQRTT(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a prior structure with the named parameters altered
%    with the specified values.
%
%    Parameters for Student-t prior [default]
%      mu       - location [0]
%      s2       - scale [1]
%      nu       - degrees of freedom [4]
%      mu_prior - prior for mu [prior_fixed]
%      s2_prior - prior for s2 [prior_fixed]
%      nu_prior - prior for nu [prior_fixed]
%
%  See also
%    PRIOR_t, PRIOR_*

% Copyright (c) 2000-2001,2010 Aki Vehtari
% Copyright (c) 2009 Jarno Vanhatalo
% Copyright (c) 2010 Jaakko Riihim�ki

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_SQRTT';
  ip.addOptional('p', [], @isstruct);
  ip.addParamValue('mu',0, @(x) isscalar(x));
  ip.addParamValue('mu_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('s2',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('s2_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('nu',4, @(x) isscalar(x) && x>0);
  ip.addParamValue('nu_prior',[], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Sqrt-t';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Sqrt-t')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('mu',ip.UsingDefaults)
    p.mu = ip.Results.mu;
  end
  if init || ~ismember('s2',ip.UsingDefaults)
    p.s2 = ip.Results.s2;
  end
  if init || ~ismember('nu',ip.UsingDefaults)
    p.nu = ip.Results.nu;
  end
  % Initialize prior structure
  if init
    p.p=[];
  end
  if init || ~ismember('mu_prior',ip.UsingDefaults)
    p.p.mu=ip.Results.mu_prior;
  end
  if init || ~ismember('s2_prior',ip.UsingDefaults)
    p.p.s2=ip.Results.s2_prior;
  end
  if init || ~ismember('nu_prior',ip.UsingDefaults)
    p.p.nu=ip.Results.nu_prior;
  end

  if init
    % set functions
    p.fh.pak = @prior_sqrtt_pak;
    p.fh.unpak = @prior_sqrtt_unpak;
    p.fh.lp = @prior_sqrtt_lp;
    p.fh.lpg = @prior_sqrtt_lpg;
    p.fh.recappend = @prior_sqrtt_recappend;
  end

end

function [w, s] = prior_sqrtt_pak(p)
  
  w=[];
  s={};
  if ~isempty(p.p.mu)
    w = p.mu;
    s=[s; 'Sqrt-Student-t.mu'];
  end        
  if ~isempty(p.p.s2)
    w = [w log(p.s2)];
    s=[s; 'log(Sqrt-Student-t.s2)'];
  end
  if ~isempty(p.p.nu)
    w = [w log(p.nu)];
    s=[s; 'log(Sqrt-Student-t.nu)'];
  end
end

function [p, w] = prior_sqrtt_unpak(p, w)
  
  if ~isempty(p.p.mu)
    i1=1;
    p.mu = w(i1);
    w = w(i1+1:end);
  end
  if ~isempty(p.p.s2)
    i1=1;
    p.s2 = exp(w(i1));
    w = w(i1+1:end);
  end
  if ~isempty(p.p.nu)
    i1=1;
    p.nu = exp(w(i1));
    w = w(i1+1:end);
  end
end

function lp = prior_sqrtt_lp(x, p)
  
  lJ = -log(2*sqrt(x));  % log(1/(2*x^(1/2))) log(|J|) of transformation
  xt = sqrt(x);          % transformation
  lp = sum(gammaln((p.nu+1)./2) - gammaln(p.nu./2) - 0.5*log(p.nu.*pi.*p.s2) - (p.nu+1)./2.*log(1+(xt-p.mu).^2./p.nu./p.s2) + lJ);
  
  if ~isempty(p.p.mu)
    lp = lp + p.p.mu.fh.lp(p.mu, p.p.mu);
  end
  if ~isempty(p.p.s2)
    lp = lp + p.p.s2.fh.lp(p.s2, p.p.s2) + log(p.s2);
  end
  if ~isempty(p.p.nu)
    lp = lp + p.p.nu.fh.lp(p.nu, p.p.nu) + log(p.nu);
  end
end

function lpg = prior_sqrtt_lpg(x, p)

  lJg = -1./(2*x);        % gradient of log(|J|) of transformation
  xt  = sqrt(x);          % transformation
  xtg = 1./(2*sqrt(x));   % derivative of transformation
  lpg = xtg.*(-(p.nu+1).*(xt-p.mu)./(p.nu.*p.s2 + (xt-p.mu).^2)) + lJg;
  
  if ~isempty(p.p.mu)
    lpgmu = sum((p.nu+1).*(xt-p.mu)./(p.nu.*p.s2 + (xt-p.mu).^2)) + p.p.mu.fh.lpg(p.mu, p.p.mu);
    lpg = [lpg lpgmu];
  end
  if ~isempty(p.p.s2)
    lpgs2 = (sum(-1./(2.*p.s2)+((p.nu + 1)*(p.mu - xt)^2)./(2*p.s2*((p.mu-xt)^2 + p.nu*p.s2))) + p.p.s2.fh.lpg(p.s2, p.p.s2)).*p.s2 + 1;
    lpg = [lpg lpgs2];
  end
  if ~isempty(p.p.nu)
    lpgnu = (0.5*sum(digamma1((p.nu+1)./2)-digamma1(p.nu./2)-1./p.nu-log(1+(xt-p.mu).^2./p.nu./p.s2)+(p.nu+1)./(1+(xt-p.mu).^2./p.nu./p.s2).*(xt-p.mu).^2./p.s2./p.nu.^2) + p.p.nu.fh.lpg(p.nu, p.p.nu)).*p.nu + 1;
    lpg = [lpg lpgnu];
  end
end

function rec = prior_sqrtt_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
  if ~isempty(p.p.mu)
    rec.mu(ri,:) = p.mu;
  end        
  if ~isempty(p.p.s2)
    rec.s2(ri,:) = p.s2;
  end
  if ~isempty(p.p.nu)
    rec.nu(ri,:) = p.nu;
  end
end
