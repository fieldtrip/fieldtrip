function p = prior_loggaussian(varargin)
%PRIOR_LOGGAUSSIAN  Log-Gaussian prior structure     
%       
%  Description
%    P = PRIOR_LOGGAUSSIAN('PARAM1', VALUE1, 'PARAM2', VALUE2, ...) 
%    creates Log-Gaussian prior structure in which the named
%    parameters have the specified values. Any unspecified
%    parameters are set to default values.
%    
%    P = PRIOR_LOGGAUSSIAN(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a prior structure with the named parameters altered
%    with the specified values.
%
%    Parameters for Log-Gaussian prior [default]
%      mu       - location [0]
%      s2       - scale squared (variance) [1]
%      mu_prior - prior for mu [prior_fixed]
%      s2_prior - prior for s2 [prior_fixed]
%
%  See also
%    PRIOR_*

% Copyright (c) 2000-2001,2010 Aki Vehtari
% Copyright (c) 2010 Jaakko Riihim�ki

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_LOGGAUSSIAN';
  ip.addOptional('p', [], @isstruct);
  ip.addParamValue('mu',0, @(x) isscalar(x));
  ip.addParamValue('mu_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('s2',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('s2_prior',[], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Log-Gaussian';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Log-Gaussian')
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

  if init
    % set functions
    p.fh.pak = @prior_loggaussian_pak;
    p.fh.unpak = @prior_loggaussian_unpak;
    p.fh.lp = @prior_loggaussian_lp;
    p.fh.lpg = @prior_loggaussian_lpg;
    p.fh.recappend = @prior_loggaussian_recappend;
  end

end

function [w, s] = prior_loggaussian_pak(p)
  
  w=[];
  s={};
  if ~isempty(p.p.mu)
    w = p.mu;
    s=[s; 'Log-Gaussian.mu'];
  end
  if ~isempty(p.p.s2)
    w = [w log(p.s2)];
    s=[s; 'log(Log-Gaussian.s2)'];
  end
end

function [p, w] = prior_loggaussian_unpak(p, w)

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
end

function lp = prior_loggaussian_lp(x, p)
  
  lJ = -log(x);    % =log(1/x)=log(|J|) of transformation
  xt  = log(x);    % transformed x
  lp = 0.5*sum(-log(2*pi) -log(p.s2)- 1./p.s2 .* sum((xt-p.mu).^2,1)) +lJ;
    
  if ~isempty(p.p.mu)
    lp = lp + p.p.mu.fh.lp(p.mu, p.p.mu);
  end
  if ~isempty(p.p.s2)
    lp = lp + p.p.s2.fh.lp(p.s2, p.p.s2) + log(p.s2);
  end
end

function lpg = prior_loggaussian_lpg(x, p)
  
  lJg = -1./x;     % gradient of log(|J|) of transformation
  xt   = log(x);   % transformed x
  xtg  = 1./x;     % derivative of transformation
  lpg = xtg.*(1./p.s2).*(p.mu-xt) + lJg;
  
  if ~isempty(p.p.mu)
    lpgmu = sum((1./p.s2).*(xt-p.mu)) + p.p.mu.fh.lpg(p.mu, p.p.mu);
    lpg = [lpg lpgmu];
  end
  if ~isempty(p.p.s2)
    lpgs2 = (sum(-0.5*(1./p.s2-1./p.s2.^2.*(xt-p.mu).^2 )) + p.p.s2.fh.lpg(p.s2, p.p.s2)).*p.s2 + 1;
    lpg = [lpg lpgs2];
  end
end

function rec = prior_loggaussian_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
  if ~isempty(p.p.mu)
    rec.mu(ri,:) = p.mu;
  end
  if ~isempty(p.p.s2)
    rec.s2(ri,:) = p.s2;
  end
end    
