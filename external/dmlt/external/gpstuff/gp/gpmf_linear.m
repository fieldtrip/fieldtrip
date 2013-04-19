function gpmf = gpmf_linear(varargin)
%GPMF_LINEAR  Create a linear mean function
%
%  Description
%    GPMF = GPMF_LINEAR('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates linear mean function structure in which the named
%    parameters have the specified values. Any unspecified
%    parameters are set to default values.
%
%    GPMF = GPMF_LINEAR(GPMF,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a mean function structure with the named parameters
%    altered with the specified values.
%  
%    Parameters for linear mean function
%      prior_mean        - prior mean (scalar or vector) for base
%                          functions' weight prior (default 0)
%      prior_cov         - prior covariances (scalar or vector) 
%                          for base functions' prior corresponding
%                          each selected input dimension. In 
%                          multiple dimension case prior_cov is a
%                          struct containing scalars or vectors.
%                          The covariances must all be either
%                          scalars (diagonal cov.matrix) or
%                          vectors (for non-diagonal cov.matrix)
%                          (default 100)  
%      selectedVariables - vector defining which inputs are active
% 
%  See also
%    GP_SET, GPMF_CONSTANT, GPMF_SQUARED, 

% Copyright (c) 2010 Tuomas Nikoskinen
% Copyright (c) 2011 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPMF_LINEAR';
  ip.addOptional('gpmf', [], @isstruct);
  ip.addParamValue('selectedVariables',[], @(x) isvector(x) && all(x>0));
  ip.addParamValue('prior_mean',0, @(x) isvector(x));
  ip.addParamValue('prior_cov',100, @(x) isvector(x));
  ip.addParamValue('mean_prior', [], @isstruct);
  ip.addParamValue('cov_prior', [], @isstruct);
  ip.parse(varargin{:});
  gpmf=ip.Results.gpmf;
  
  if isempty(gpmf)
    % Initialize a mean function
    init=true;
    gpmf.type = 'gpmf_linear';
  else
    % Modify a mean function
    if ~isfield(gpmf,'type') && isequal(gpmf.type,'gpmf_linear')
      error('First argument does not seem to be a linear mean function')
    end
    init=false;
  end
  % Initialize parameters
  if init || ~ismember('prior_mean',ip.UsingDefaults)
    gpmf.b=ip.Results.prior_mean(:)';
  end
  if init || ~ismember('prior_cov',ip.UsingDefaults)
    gpmf.B=ip.Results.prior_cov(:)';
  end
  if ~ismember('selectedVariables',ip.UsingDefaults)
    gpmf.selectedVariables=ip.Results.selectedVariables;
  end
  if init || ~ismember('mean_prior',ip.UsingDefaults)
    gpmf.p.b=ip.Results.mean_prior;
  end
  if init || ~ismember('cov_prior',ip.UsingDefaults)
    gpmf.p.B=ip.Results.cov_prior;
  end
  if init
    % Set the function handles to the nested functions
    gpmf.fh.geth = @gpmf_geth;
    gpmf.fh.pak = @gpmf_pak;
    gpmf.fh.unpak = @gpmf_unpak;
    gpmf.fh.lp = @gpmf_lp;
    gpmf.fh.lpg = @gpmf_lpg;
    gpmf.fh.recappend = @gpmf_recappend;
  end

end

function h = gpmf_geth(gpmf, x)
%GPMF_GETH  Calculate the base function values for given input.
%
%  Description
%    H = GPMF_GETH(GPMF,X) takes in a mean function structure
%    GPMF and inputs X. The function returns the linear base
%    function values H in the given input points. If
%    selectedVariables is used the function returns only the
%    values corresponding active inputs. The base function values
%    are returned as a matrix in which each row corresponds to
%    one dimension and the first row is for the smallest
%    dimension.
  
  if ~isfield(gpmf,'selectedVariables')
    h=x';
  else
    h=x(:,gpmf.selectedVariables)';
  end
  
end

function [w, s] = gpmf_pak(gpmf, w)
%GPMF_PAK  Combine GP mean function parameters into one vector
%
%  Description
%
%  See also

  
  w = []; s = {};
  if ~isempty(gpmf.p.b)
    w = gpmf.b;
    if numel(gpmf.b)>1
      s = [s; sprintf('gpmf_linear.b x %d',numel(gpmf.b))];
    else
      s = [s; 'gpmf_linear.b'];
    end
    % Hyperparameters of coeffSigma2
    [wh sh] = gpmf.p.b.fh.pak(gpmf.p.b);
    w = [w wh];
    s = [s; sh];
  end
  
  if ~isempty(gpmf.p.B)
    w = [w log(gpmf.B)];
    if numel(gpmf.B)>1
      s = [s; sprintf('log(gpmf_linear.B x %d)',numel(gpmf.B))];
    else
      s = [s; 'log(gpmf_linear.B)'];
    end
    % Hyperparameters of coeffSigma2
    [wh sh] = gpmf.p.B.fh.pak(gpmf.p.B);
    w = [w wh];
    s = [s; sh];
  end
  
end

function [gpmf, w] = gpmf_unpak(gpmf, w)
%GPCF_LINEAR_UNPAK  Sets the mean function parameters 
%                   into the structure
%
%  Description
%
%  See also
  
  gpp=gpmf.p;

  if ~isempty(gpp.b)
    i2=length(gpmf.b);
    i1=1;
    gpmf.b = w(i1:i2);
    w = w(i2+1:end);
    
    % Hyperparameters of b
    [p, w] = gpmf.p.b.fh.unpak(gpmf.p.b, w);
    gpmf.p.b = p;
  end
  
  if ~isempty(gpp.B)
    i2=length(gpmf.B);
    i1=1;
    gpmf.B = exp(w(i1:i2));
    w = w(i2+1:end);
    
    % Hyperparameters of b
    [p, w] = gpmf.p.B.fh.unpak(gpmf.p.B, w);
    gpmf.p.B = p;
  end
  
end

function lp = gpmf_lp(gpmf)
%GPCF_SEXP_LP  Evaluate the log prior of covariance function parameters
%
%  Description
%
%  See also

% Evaluate the prior contribution to the error. The parameters that
% are sampled are transformed, e.g., W = log(w) where w is all
% the "real" samples. On the other hand errors are evaluated in
% the W-space so we need take into account also the Jacobian of
% transformation, e.g., W -> w = exp(W). See Gelman et.al., 2004,
% Bayesian data Analysis, second edition, p24.
  lp = 0;
  gpp=gpmf.p;
  
  if ~isempty(gpmf.p.b)
    lp = lp + gpp.b.fh.lp(gpmf.b, ...
                    gpp.b);
  end

  if ~isempty(gpp.B)
    lp = lp + gpp.B.fh.lp(gpmf.B, ...
                    gpp.B) +sum(log(gpmf.B));
  end
end

function [lpg_b, lpg_B] = gpmf_lpg(gpmf)
%GPCF_SEXP_LPG  Evaluate gradient of the log prior with respect
%               to the parameters.
%
%  Description
%    LPG = GPCF_SEXP_LPG(GPCF) takes a covariance function
%    structure GPCF and returns LPG = d log (p(th))/dth, where th
%    is the vector of parameters.
%
%  See also
%    GPCF_SEXP_PAK, GPCF_SEXP_UNPAK, GPCF_SEXP_LP, GP_G

  lpg_b=[];, lpg_B=[];
  gpp=gpmf.p;
  
  if ~isempty(gpmf.p.b)
    lll = length(gpmf.b);
    lpgs = gpp.b.fh.lpg(gpmf.b, gpp.b);
    lpg_b = [lpgs(1:lll) lpgs(lll+1:end)]; %.*gpmf.b+1
  end
  
  if ~isempty(gpmf.p.B)
    lll = length(gpmf.B);
    lpgs = gpp.B.fh.lpg(gpmf.B, gpp.B);
    lpg_B = [lpgs(1:lll).*gpmf.B+1 lpgs(lll+1:end)];
  end
end

function recmf = gpmf_recappend(recmf, ri, gpmf)
%RECAPPEND  Record append
%
%  Description
%
%  See also
%    GP_MC and GP_MC -> RECAPPEND

% Initialize record
  if nargin == 2
    recmf.type = 'gpmf_linear';

    % Initialize parameters
    recmf.b= [];
    recmf.B = [];

    % Set the function handles
    recmf.fh.geth = @gpmf_geth;
    recmf.fh.pak = @gpmf_pak;
    recmf.fh.unpak = @gpmf_unpak;
    recmf.fh.lp = @gpmf_lp;
    recmf.fh.lpg = @gpmf_lpg;
    recmf.fh.recappend = @gpmf_recappend;

    recmf.p=[];
    recmf.p.b=[];
    recmf.p.B=[];
    if isfield(ri.p,'b') && ~isempty(ri.p.b)
      recmf.p.b = ri.p.b;
    end
    if ~isempty(ri.p.B)
      recmf.p.B = ri.p.B;
    end
    return
  end

  gpp = gpmf.p;

  % record magnSigma2
  if ~isempty(gpmf.b)
    recmf.b(ri,:)=gpmf.b;
    if ~isempty(recmf.p.b)
      recmf.p.b = gpp.b.fh.recappend(recmf.p.b, ri, gpmf.p.b);
    end
  elseif ri==1
    recmf.b=[];
  end
  
  if ~isempty(gpmf.B)
    recmf.B(ri,:)=gpmf.B;
    if ~isempty(recmf.p.B)
      recmf.p.B = gpp.B.fh.recappend(recmf.p.B, ri, gpmf.p.B);
    end
  elseif ri==1
    recmf.B=[];
  end

end
