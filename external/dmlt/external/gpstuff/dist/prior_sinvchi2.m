function p = prior_sinvchi2(varargin)
%PRIOR_SINVCHI2  Scaled-Inv-Chi^2 prior structure
%       
%  Description
%    P = PRIOR_SINVCHI2('PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    creates Scaled-Inv-Chi^2 prior structure in which the named
%    parameters have the specified values. Any unspecified
%    parameters are set to default values.
%
%    P = PRIOR_SINVCHI2(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a prior structure with the named parameters altered
%    with the specified values.
%
%    Parameterisation is done by Bayesian Data Analysis,  
%    second edition, Gelman et.al 2004.
%    
%    Parameters for Scaled-Inv-Chi^2 [default]
%      s2       - scale squared (variance) [1]
%      nu       - degrees of freedom [4]
%      s2_prior - prior for s2 [prior_fixed]
%      nu_prior - prior for nu [prior_fixed]
%  
%  See also
%    PRIOR_*

% Copyright (c) 2000-2001,2010 Aki Vehtari
% Copyright (c) 2010 Jaakko Riihim�ki

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_SINVCHI2';
  ip.addOptional('p', [], @isstruct);
  ip.addParamValue('s2',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('s2_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('nu',4, @(x) isscalar(x) && x>0);
  ip.addParamValue('nu_prior',[], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'S-Inv-Chi2';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'S-Inv-Chi2')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end

  % Initialize parameters
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
  if init || ~ismember('s2_prior',ip.UsingDefaults)
    p.p.s2=ip.Results.s2_prior;
  end
  if init || ~ismember('nu_prior',ip.UsingDefaults)
    p.p.nu=ip.Results.nu_prior;
  end

  if init
    % set functions
    p.fh.pak = @prior_sinvchi2_pak;
    p.fh.unpak = @prior_sinvchi2_unpak;
    p.fh.lp = @prior_sinvchi2_lp;
    p.fh.lpg = @prior_sinvchi2_lpg;
    p.fh.recappend = @prior_sinvchi2_recappend;
  end

end

function [w, s] = prior_sinvchi2_pak(p)
  
  w=[];
  s={};
  if ~isempty(p.p.s2)
    w = log(p.s2);
    s=[s; 'log(Sinvchi2.s2)'];
  end
  if ~isempty(p.p.nu)
    w = [w log(p.nu)];
    s=[s; 'log(Sinvchi2.nu)'];
  end
end

function [p, w] = prior_sinvchi2_unpak(p, w)
  
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

function lp = prior_sinvchi2_lp(x, p)
  lp = -sum((p.nu./2+1) .* log(x) + (p.s2.*p.nu./2./x) + (p.nu/2) .* log(2./(p.s2.*p.nu)) + gammaln(p.nu/2)) ;
  
  if ~isempty(p.p.s2)
    lp = lp + p.p.s2.fh.lp(p.s2, p.p.s2) + log(p.s2);
  end
  if ~isempty(p.p.nu)
    lp = lp + p.p.nu.fh.lp(p.nu, p.p.nu) + log(p.nu);
  end
end

function lpg = prior_sinvchi2_lpg(x, p)
  lpg = -(p.nu/2+1)./x +p.nu.*p.s2./(2*x.^2);

  if ~isempty(p.p.s2)
    lpgs2 = (-sum(p.nu/2.*(1./x-1./p.s2)) + p.p.s2.fh.lpg(p.s2, p.p.s2)).*p.s2 + 1; 
    lpg = [lpg lpgs2];
  end
  if ~isempty(p.p.nu)
    lpgnu = (-sum(0.5*(log(x) + p.s2./x + log(2./p.s2./p.nu) - 1 + digamma1(p.nu/2))) + p.p.nu.fh.lpg(p.nu, p.p.nu)).*p.nu + 1;
    lpg = [lpg lpgnu];
  end
end

function rec = prior_sinvchi2_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
  if ~isempty(p.p.s2)
    rec.s2(ri,:) = p.s2;
  end
  if ~isempty(p.p.nu)
    rec.nu(ri,:) = p.nu;
  end
end
