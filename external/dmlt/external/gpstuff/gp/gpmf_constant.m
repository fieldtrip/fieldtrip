function gpmf = gpmf_constant(varargin)
%GPMF_CONSTANT  Create a constant mean function
%
%  Description
%    GPMF = GPMF_CONSTANT('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates constant mean function structure in which the named
%    parameters have the specified values. Any unspecified
%    parameters are set to default values.
%
%    GPMF = GPMF_CONSTANT(GPMF,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a mean function structure with the named parameters
%    altered with the specified values.
%  
%    Parameters for constant mean function
%      constant          - constant value for the constant
%                          base function (default 1)
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
% 
%  See also
%    GP_SET, GPMF_LINEAR, GPMF_SQUARED
%
  
% Copyright (c) 2010 Tuomas Nikoskinen
% Copyright (c) 2011 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPMF_CONSTANT';
  ip.addOptional('gpmf', [], @isstruct);
  ip.addParamValue('constant',1, @(x) isvector(x) && all(x>0));
  ip.addParamValue('prior_mean',0, @(x) isvector(x));
  ip.addParamValue('prior_cov',100, @(x) isvector(x));
  ip.addParamValue('mean_prior', [], @isstruct);
  ip.addParamValue('cov_prior', [], @isstruct);
  ip.parse(varargin{:});
  gpmf=ip.Results.gpmf;
  
  if isempty(gpmf)
    % Initialize a mean function
    init=true;
    gpmf.type = 'gpmf_constant';
  else
    % Modify a mean function
    if ~isfield(gpmf,'type') && isequal(gpmf.type,'gpmf_constant')
      error('First argument does not seem to be a constant mean function')
    end
    init=false;
  end
  
  % Initialize parameters
  if init || ~ismember('type',ip.UsingDefaults)
    gpmf.constant = ip.Results.constant;
  end
  if init || ~ismember('prior_mean',ip.UsingDefaults)
    gpmf.b=ip.Results.prior_mean(:)';
  end
  if init || ~ismember('prior_mean',ip.UsingDefaults)
    gpmf.B=ip.Results.prior_cov(:)';
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
%GPMF_GETH  Calculate the base function values for a given input.
%
%  Description
%    H = GPMF_GETH(GPMF,X) takes in a mean function structure
%    GPMF and inputs X. The function returns a row vector of
%    length(X) containing the constant value which is by default
%    1.
  
  constant=gpmf.constant;
  h = repmat(constant,1,length(x(:,1)));
  
end

function [w, s] = gpmf_pak(gpmf, w)
%GPMF_PAK  Combine GP mean function parameters into one vector
%
%  Description
%    W = GPCF_LINEAR_PAK(GPCF) takes a covariance function
%    structure GPCF and combines the covariance function
%    parameters and their hyperparameters into a single row
%    vector W.
%
%       w = [ log(gpcf.coeffSigma2)
%             (hyperparameters of gpcf.coeffSigma2)]'
%
%  See also
%    GPCF_LINEAR_UNPAK
  
  w = []; s = {};
  if ~isempty(gpmf.p.b)
    w = gpmf.b;
    if numel(gpmf.b)>1
      s = [s; sprintf('gpmf_constant.b x %d',numel(gpmf.b))];
    else
      s = [s; 'gpmf_constant.b'];
    end
    % Hyperparameters of coeffSigma2
    [wh sh] = gpmf.p.b.fh.pak(gpmf.p.b);
    w = [w wh];
    s = [s; sh];
  end
  
  if ~isempty(gpmf.p.B)
    w = [w log(gpmf.B)];
    if numel(gpmf.B)>1
      s = [s; sprintf('log(gpmf_constant.B x %d)',numel(gpmf.B))];
    else
      s = [s; 'log(gpmf_constant.B)'];
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
%    [GPCF, W] = GPMF_UNPAK(GPCF, W) takes a covariance
%    function structure GPCF and a hyper-parameter vector W, and
%    returns a covariance function structure identical to the
%    input, except that the covariance hyper-parameters have been
%    set to the values in W. Deletes the values set to GPCF from
%    W and returns the modified W.
%
%    Assignment is inverse of  
%       w = [ log(gpcf.coeffSigma2)
%             (hyperparameters of gpcf.coeffSigma2)]'
%
%  See also
%   GPCF_LINEAR_PAK
  
  gpp=gpmf.p;

  if ~isempty(gpp.b)
    i2=length(gpmf.b);
    i1=1;
    gpmf.b = w(i1:i2);
    w = w(i2+1:end);
    
    % Hyperparameters of coeffSigma2
    [p, w] = gpmf.p.b.fh.unpak(gpmf.p.b, w);
    gpmf.p.b = p;
  end
  
  if ~isempty(gpp.B)
    i2=length(gpmf.B);
    i1=1;
    gpmf.B = exp(w(i1:i2));
    w = w(i2+1:end);
    
    % Hyperparameters of coeffSigma2
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

  lpg_b=[]; lpg_B=[];
  gpp=gpmf.p;
  
  if ~isempty(gpmf.p.b)
    lll = length(gpmf.b);
    lpgs = gpp.b.fh.lpg(gpmf.b, gpp.b);
    lpg_b = [lpgs(1:lll) lpgs(lll+1:end)];  %.*gpmf.b+1
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
    recmf.type = 'gpmf_constant';

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
