function p = prior_invgamma(varargin)
%PRIOR_INVGAMMA  Inverse-gamma prior structure     
%       
%  Description
%    P = PRIOR_INVGAMMA('PARAM1', VALUE1, 'PARAM2', VALUE2, ...) 
%    creates Gamma prior structure in which the named parameters
%    have the specified values. Any unspecified parameters are set
%    to default values.
%
%    P = PRIOR_INVGAMMA(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a prior structure with the named parameters altered
%    with the specified values.
%  
%    Parameterisation is done by Bayesian Data Analysis,  
%    second edition, Gelman et.al 2004.
%
%    Parameters for Gamma prior [default]
%      sh       - shape [4]
%      s        - scale [1]
%      sh_prior - prior for sh [prior_fixed]
%      s_prior  - prior for s [prior_fixed]
%
%  See also
%    PRIOR_*


% Copyright (c) 2000-2001,2010 Aki Vehtari
% Copyright (c) 2010 Jaakko Riihim�ki

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_INVGAMMA';
  ip.addOptional('p', [], @isstruct);
  ip.addParamValue('sh',4, @(x) isscalar(x) && x>0);
  ip.addParamValue('sh_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('s',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('s_prior',[], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Inv-Gamma';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Inv-Gamma')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('sh',ip.UsingDefaults)
    p.sh = ip.Results.sh;
  end
  if init || ~ismember('s',ip.UsingDefaults)
    p.s = ip.Results.s;
  end
  % Initialize prior structure
  if init
    p.p=[];
  end
  if init || ~ismember('sh_prior',ip.UsingDefaults)
    p.p.sh=ip.Results.sh_prior;
  end
  if init || ~ismember('s_prior',ip.UsingDefaults)
    p.p.s=ip.Results.s_prior;
  end

  if init
    % set functions
    p.fh.pak = @prior_invgamma_pak;
    p.fh.unpak = @prior_invgamma_unpak;
    p.fh.lp = @prior_invgamma_lp;
    p.fh.lpg = @prior_invgamma_lpg;
    p.fh.recappend = @prior_invgamma_recappend;
  end

end

function [w, s] = prior_invgamma_pak(p)
  
  w=[];
  s={};
  if ~isempty(p.p.sh)
    w = log(p.sh);
    s=[s; 'log(Invgamma.sh)'];
  end
  if ~isempty(p.p.s)
    w = [w log(p.s)];
    s=[s; 'log(Invgamma.s)'];
  end
end

function [p, w] = prior_invgamma_unpak(p, w)

  if ~isempty(p.p.sh)
    i1=1;
    p.sh = exp(w(i1));
    w = w(i1+1:end);
  end
  if ~isempty(p.p.s)
    i1=1;
    p.s = exp(w(i1));
    w = w(i1+1:end);
  end
end

function lp = prior_invgamma_lp(x, p)
  
  lp = sum(-p.s./x - (p.sh+1).*log(x) +p.sh.*log(p.s)  - gammaln(p.sh));
  
  if ~isempty(p.p.sh)
    lp = lp + p.p.sh.fh.lp(p.sh, p.p.sh) + log(p.sh);
  end
  if ~isempty(p.p.s)
    lp = lp + p.p.s.fh.lp(p.s, p.p.s) + log(p.s);
  end
end

function lpg = prior_invgamma_lpg(x, p)
  
  lpg = -(p.sh+1)./x + p.s./x.^2;
  
  if ~isempty(p.p.sh)
    lpgsh = (-sum(digamma1(p.sh) + log(p.s) - log(x) ) + p.p.sh.fh.lpg(p.sh, p.p.sh)).*p.sh + 1;
    lpg = [lpg lpgsh];
  end
  if ~isempty(p.p.s)
    lpgs = (sum(p.sh./p.s+1./x) + p.p.s.fh.lpg(p.s, p.p.s)).*p.s + 1;
    lpg = [lpg lpgs];
  end
  
end

function rec = prior_invgamma_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
  if ~isempty(p.p.sh)
    rec.sh(ri,:) = p.sh;
  end
  if ~isempty(p.p.s)
    rec.s(ri,:) = p.s;
  end
end
