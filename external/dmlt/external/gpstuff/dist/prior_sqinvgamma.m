function p = prior_sqinvgamma(varargin)
%PRIOR_SQINVGAMMA  Gamma prior structure for square inverse of the parameter
%
%  Description
%    P = PRIOR_SQINVGAMMA('PARAM1', VALUE1, 'PARAM2', VALUE2, ...) 
%    creates Gamma prior structure for square inverse of the
%    parameter in which the named parameters have the specified
%    values. Any unspecified parameters are set to default values.
%
%    P = PRIOR_SQINVGAMMA(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a prior structure with the named parameters altered
%    with the specified values.
%  
%    Parametrisation is done by Bayesian Data Analysis,  
%    second edition, Gelman et.al. 2004.
%
%    Parameters for Gamma prior [default]
%      sh       - shape [4]
%      is       - inverse scale [1]
%      sh_prior - prior for sh [prior_fixed]
%      is_prior - prior for is [prior_fixed]
%
%  See also
%    PRIOR_*

% Copyright (c) 2000-2001,2010,2012 Aki Vehtari
% Copyright (c) 2010 Jaakko Riihim�ki

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_SQINVGAMMA';
  ip.addOptional('p', [], @isstruct);
  ip.addParamValue('sh',4, @(x) isscalar(x) && x>0);
  ip.addParamValue('sh_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('is',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('is_prior',[], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'SqInv-Gamma';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'SqInv-Gamma')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('sh',ip.UsingDefaults)
    p.sh = ip.Results.sh;
  end
  if init || ~ismember('is',ip.UsingDefaults)
    p.is = ip.Results.is;
  end
  % Initialize prior structure
  if init
    p.p=[];
  end
  if init || ~ismember('sh_prior',ip.UsingDefaults)
    p.p.sh=ip.Results.sh_prior;
  end
  if init || ~ismember('is_prior',ip.UsingDefaults)
    p.p.is=ip.Results.is_prior;
  end

  if init
    % set functions
    p.fh.pak = @prior_sqinvgamma_pak;
    p.fh.unpak = @prior_sqinvgamma_unpak;
    p.fh.lp = @prior_sqinvgamma_lp;
    p.fh.lpg = @prior_sqinvgamma_lpg;
    p.fh.recappend = @prior_sqinvgamma_recappend;
  end

end

function [w, s] = prior_sqinvgamma_pak(p)
  
  w=[];
  s={};
  if ~isempty(p.p.sh)
    w = log(p.sh);
    s=[s; 'log(SqInv-Gamma.sh)'];
  end
  if ~isempty(p.p.is)
    w = [w log(p.is)];
    s=[s; 'log(SqInv-Gamma.is)'];
  end
end

function [p, w] = prior_sqinvgamma_unpak(p, w)

  if ~isempty(p.p.sh)
    i1=1;
    p.sh = exp(w(i1));
    w = w(i1+1:end);
  end
  if ~isempty(p.p.is)
    i1=1;
    p.is = exp(w(i1));
    w = w(i1+1:end);
  end
end

function lp = prior_sqinvgamma_lp(x, p)
  
  lJ = -log(x)*3 + log(2);  % log(-2/x^3) log(|J|) of transformation
  xt = x.^-2;               % transformation
  lp = sum(-p.is.*xt + (p.sh-1).*log(xt) +p.sh.*log(p.is)  -gammaln(p.sh) +lJ);
  
  if ~isempty(p.p.sh)
    lp = lp + p.p.sh.fh.lp(p.sh, p.p.sh) + log(p.sh);
  end
  if ~isempty(p.p.is)
    lp = lp + p.p.is.fh.lp(p.is, p.p.is) + log(p.is);
  end
end

function lpg = prior_sqinvgamma_lpg(x, p)
  
  lJg = -3./x;              % gradient of log(|J|) of transformation
  xt  = x.^-2;              % transformation
  xtg = -2/x.^3;            % derivative of transformation
  lpg = xtg.*((p.sh-1)./xt - p.is) + lJg;
  
  if ~isempty(p.p.sh)
    lpgsh = (sum(-digamma1(p.sh) + log(p.is) + log(x)) + p.p.sh.fh.lpg(p.sh, p.p.sh)).*p.sh + 1;
    lpg = [lpg lpgsh];
  end
  if ~isempty(p.p.is)
    lpgis = (sum(p.sh./p.is+x) + p.p.is.fh.lpg(p.is, p.p.is)).*p.is + 1;
    lpg = [lpg lpgis];
  end
  
end

function rec = prior_sqinvgamma_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
  if ~isempty(p.p.sh)
    rec.sh(ri,:) = p.sh;
  end
  if ~isempty(p.p.is)
    rec.is(ri,:) = p.is;
  end
end    
