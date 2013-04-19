function p = prior_laplace(varargin)
%PRIOR_LAPLACE  Laplace (double exponential) prior structure     
%       
%  Description
%    P = PRIOR_LAPLACE('PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    creates Laplace prior structure in which the named parameters
%    have the specified values. Any unspecified parameters are set
%    to default values.
%    
%    P = PRIOR_LAPLACE(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a prior structure with the named parameters altered
%    with the specified values.
%
%    Parameters for Laplace prior [default]
%      mu       - location [0]
%      s        - scale [1]
%      mu_prior - prior for mu [prior_fixed]
%      s_prior  - prior for s  [prior_fixed]
%  
%  See also
%    PRIOR_*

% Copyright (c) 2000-2001,2010 Aki Vehtari
% Copyright (c) 2010 Jaakko Riihim�ki

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_LAPLACE';
  ip.addOptional('p', [], @isstruct);
  ip.addParamValue('mu',0, @(x) isscalar(x) && x>0);
  ip.addParamValue('mu_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('s',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('s_prior',[], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Laplace';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Laplace')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('mu',ip.UsingDefaults)
    p.mu = ip.Results.mu;
  end
  if init || ~ismember('s',ip.UsingDefaults)
    p.s = ip.Results.s;
  end
  % Initialize prior structure
  if init
    p.p=[];
  end
  if init || ~ismember('mu_prior',ip.UsingDefaults)
    p.p.mu=ip.Results.mu_prior;
  end
  if init || ~ismember('s_prior',ip.UsingDefaults)
    p.p.s=ip.Results.s_prior;
  end

  if init
    % set functions
    p.fh.pak = @prior_laplace_pak;
    p.fh.unpak = @prior_laplace_unpak;
    p.fh.lp = @prior_laplace_lp;
    p.fh.lpg = @prior_laplace_lpg;
    p.fh.recappend = @prior_laplace_recappend;
  end

end

function [w, s] = prior_laplace_pak(p)
  
  w=[];
  s={};
  if ~isempty(p.p.mu)
    w = p.mu;
    s=[s; 'Laplace.mu'];
  end
  if ~isempty(p.p.s)
    w = [w log(p.s)];
    s=[s; 'log(Laplace.s)'];
  end
end

function [p, w] = prior_laplace_unpak(p, w)
  
  if ~isempty(p.p.mu)
    i1=1;
    p.mu = w(i1);
    w = w(i1+1:end);
  end
  if ~isempty(p.p.s)
    i1=1;
    p.s = exp(w(i1));
    w = w(i1+1:end);
  end
end

function lp = prior_laplace_lp(x, p)
  
  lp = sum(-log(2*p.s) - 1./p.s.* abs(x-p.mu));
  
  if ~isempty(p.p.mu)
    lp = lp + p.p.mu.fh.lp(p.mu, p.p.mu);
  end
  if ~isempty(p.p.s)
    lp = lp + p.p.s.fh.lp(p.s, p.p.s) + log(p.s);
  end
end

function lpg = prior_laplace_lpg(x, p)

  lpg = -sign(x-p.mu)./p.s; 
  
  if ~isempty(p.p.mu)
    lpgmu = sum(sign(x-p.mu)./p.s) + p.p.mu.fh.lpg(p.mu, p.p.mu);
    lpg = [lpg lpgmu];
  end
  if ~isempty(p.p.s)
    lpgs = (sum(-1./p.s +1./p.s.^2.*abs(x-p.mu)) + p.p.s.fh.lpg(p.s, p.p.s)).*p.s + 1;
    lpg = [lpg lpgs];
  end
end

function rec = prior_laplace_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
  if ~isempty(p.p.mu)
    rec.mu(ri,:) = p.mu;
  end
  if ~isempty(p.p.s)
    rec.s(ri,:) = p.s;
  end
end    
