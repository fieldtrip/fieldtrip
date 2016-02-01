function gpcf = gpcf_scaled(varargin)
%GPCF_SCALED  Create a scaled covariance function
%
%  Description
%    GPCF = GPCF_scaled('cf', {GPCF_1, GPCF_2, ...}) 
%    creates a scaled version of a covariance function as follows
%          GPCF_scaled = diag(x(:,scaler))*GPCF*diag(x(:,scaler))
%    where x is the matrix of inputs (see, e.g. gp_trcov).
%
%    Parameters for the scaled covariance function are [default]
%      cf        - covariance function to be scaled (compulsory)
%      scaler    - the input that is used for scaling [1]
%
%  See also
%    GP_SET, GPCF_*

% For more information on models leading to scaled covariance function see,
% for example: 
%
% GELFAND, KIM, SIRMANS, and BANERJEE (2003). Spatial Modeling With
% Spatially Varying Coefficient Processes. Journal of the American
% Statistical Association June 2003, Vol. 98, No. 462
%
  
% Copyright (c) 2009-2012 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 2 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPCF_SCALED';
  ip.addOptional('gpcf', [], @isstruct);
  ip.addParamValue('cf',[], @isstruct);
  ip.addParamValue('scaler',1, @(x) isscalar(x) && x>0);
  ip.parse(varargin{:});
  gpcf=ip.Results.gpcf;

  if isempty(gpcf)
    init=true;
    gpcf.type = 'gpcf_scaled';
  else
    if ~isfield(gpcf,'type') && ~isequal(gpcf.type,'gpcf_scaled')
      error('First argument does not seem to be a valid covariance function structure')
    end
    init=false;
  end
  
  if init || ~ismember('cf',ip.UsingDefaults)
    % Initialize parameters
    gpcf.cf = {};
    cfs=ip.Results.cf;
    if ~isempty(cfs)
        gpcf.cf{1} = cfs;
    else
      error('A covariance function has to be given in cf');
    end
  end
  
  if init || ~ismember('scaler',ip.UsingDefaults)
      gpcf.scaler = ip.Results.scaler;
  end
  
  if init
    % Set the function handles to the subfunctions
    gpcf.fh.pak = @gpcf_scaled_pak;
    gpcf.fh.unpak = @gpcf_scaled_unpak;
    gpcf.fh.lp = @gpcf_scaled_lp;
    gpcf.fh.lpg = @gpcf_scaled_lpg;
    gpcf.fh.cfg = @gpcf_scaled_cfg;
    gpcf.fh.ginput = @gpcf_scaled_ginput;
    gpcf.fh.cov = @gpcf_scaled_cov;
    gpcf.fh.trcov  = @gpcf_scaled_trcov;
    gpcf.fh.trvar  = @gpcf_scaled_trvar;
    gpcf.fh.recappend = @gpcf_scaled_recappend;
  end

end

function [w, s] = gpcf_scaled_pak(gpcf)
%GPCF_scaled_PAK  Combine GP covariance function parameters into one vector
%
%  Description
%    W = GPCF_scaled_PAK(GPCF, W) loops through all the covariance
%    functions and packs their parameters into one vector as
%    described in the respective functions. This is a mandatory 
%    subfunction used for example in energy and gradient computations.
%
%  See also
%    GPCF_scaled_UNPAK
  
  w = []; s = {};
  
  cf = gpcf.cf{1};
  [wi si] = feval(cf.fh.pak, cf);
  w = [w wi];
  s = [s; si];
end

function [gpcf, w] = gpcf_scaled_unpak(gpcf, w)
%GPCF_scaled_UNPAK  Sets the covariance function parameters into
%                 the structures
%
%  Description
%    [GPCF, W] = GPCF_scaled_UNPAK(GPCF, W) loops through all the
%    covariance functions and unpacks their parameters from W to
%    each covariance function structure. This is a mandatory 
%    subfunction used for example in energy and gradient computations.
%
%  See also
%    GPCF_scaled_PAK
%
 
    cf = gpcf.cf{1};
    [cf, w] = feval(cf.fh.unpak, cf, w);
    gpcf.cf{1} = cf;
 
end

function lp = gpcf_scaled_lp(gpcf)
%GPCF_scaled_LP  Evaluate the log prior of covariance function parameters
%
%  Description
%    LP = GPCF_scaled_LP(GPCF, X, T) takes a covariance function
%    structure GPCF and returns log(p(th)), where th collects the
%    parameters. This is a mandatory subfunction used for example 
%    in energy computations.
%
%  See also
%    GPCF_scaled_PAK, GPCF_scaled_UNPAK, GPCF_scaled_LPG, GP_E
  
  lp = 0;
  
  cf = gpcf.cf{1};
  lp = lp + feval(cf.fh.lp, cf);
  
end

function lpg = gpcf_scaled_lpg(gpcf)
%GPCF_scaled_LPG  Evaluate gradient of the log prior with respect
%               to the parameters.
%
%  Description
%    LPG = GPCF_scaled_LPG(GPCF) takes a covariance function
%    structure GPCF and returns LPG = d log (p(th))/dth, where th
%    is the vector of parameters. This is a mandatory subfunction 
%    used for example in gradient computations.
%
%  See also
%    GPCF_scaled_PAK, GPCF_scaled_UNPAK, GPCF_scaled_LP, GP_G
  lpg = [];
  
  % Evaluate the gradients
  cf = gpcf.cf{1};
  lpg=[lpg cf.fh.lpg(cf)];
end

function DKff = gpcf_scaled_cfg(gpcf, x, x2, mask, i1)
%GPCF_scaled_CFG  Evaluate gradient of covariance function
%               with respect to the parameters.
%
%  Description
%    DKff = GPCF_scaled_CFG(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to th (cell array with matrix elements). This is a 
%    mandatory subfunction used in gradient computations.
%
%    DKff = GPCF_scaled_CFG(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_scaled_CFG(GPCF, X, [], MASK) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the diagonal of gradients of covariance matrix
%    Kff = k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_scaled_CFG(GPCF, X, X2, [], i) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2), or k(X,X) if X2 is empty, with respect to ith 
%    hyperparameter. This subfunction is needed when using memory
%    save option in gp_set.
%
%  See also
%    GPCF_scaled_PAK, GPCF_scaled_UNPAK, GPCF_scaled_LP, GP_G

  [n, m] =size(x);
  if nargin==5
    % Use memory save option
    savememory=1;
    if i1==0
      % Return number of hyperparameters
      DKff=gpcf.cf{1}.fh.cfg(gpcf.cf{1},[],[],[],0);
      return
    end
  else
    savememory=0;
  end
  DKff = {};
  % Evaluate: DKff{1} = d Kff / d magnSigma2
  %           DKff{2} = d Kff / d lengthScale
  % NOTE! Here we have already taken into account that the parameters are transformed
  % through log() and thus dK/dlog(p) = p * dK/dp

  % evaluate the gradient for training covariance
  if nargin == 2 || (isempty(x2) && isempty(mask))
    
    scale = sparse(1:n,1:n,x(:,gpcf.scaler),n,n);
    
    % Evaluate the gradients
    DKff = {};
    cf = gpcf.cf{1};
    if ~savememory
      DK = cf.fh.cfg(cf, x);
    else
      DK = {cf.fh.cfg(cf,x,[],[],i1)};
    end
    for j = 1:length(DK)
        DKff{end+1} = scale*DK{j}*scale;
    end
    
    % Evaluate the gradient of non-symmetric covariance (e.g. K_fu)
  elseif nargin == 3 || isempty(mask)
    if size(x,2) ~= size(x2,2)
      error('gpcf_scaled -> _ghyper: The number of columns in x and x2 has to be the same. ')
    end
    scale = sparse(1:n,1:n,x(:,gpcf.scaler),n,n);
    n2 = length(x2);
    scale2 = sparse(1:n2,1:n2,x2(:,gpcf.scaler),n2,n2);
    
    % Evaluate the gradients
    DKff = {};
    cf = gpcf.cf{1};
    if ~savememory
      DK = cf.fh.cfg(cf, x, x2);
    else
      DK = {cf.fh.cfg(cf,x, x2, [], i1)};
    end
    
    for j = 1:length(DK)
        DKff{end+1} = scale*DK{j}*scale2;
    end
    
    % Evaluate: DKff{1}    = d mask(Kff,I) / d magnSigma2
    %           DKff{2...} = d mask(Kff,I) / d lengthScale
  elseif nargin == 4 || nargin == 5
      
      % Evaluate the gradients
      DKff = {};
      scale = x(:,gpcf.scaler);
      
      cf = gpcf.cf{1};
      if ~savememory
        DK = cf.fh.cfg(cf, x, [], 1);
      else
        DK = cf.fh.cfg(cf, x, [], 1, i1);
      end
      
      for j = 1:length(DK)
          DKff{end+1} = scale.*DK{j}.*scale;
      end
      
  end
  if savememory
    DKff=DKff{1};
  end
end


function DKff = gpcf_scaled_ginput(gpcf, x, x2,i1)
%GPCF_scaled_GINPUT  Evaluate gradient of covariance function with 
%                  respect to x
%
%  Description
%    DKff = GPCF_scaled_GINPUT(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to X (cell array with matrix elements). This subfunction 
%    is needed when computing gradients with respect to inducing 
%    inputs in sparse approximations.
%
%    DKff = GPCF_scaled_GINPUT(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_scaled_GINPUT(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2), or k(X,X) if X2 is empty, with respect to ith 
%    covariate in X (cell array with matrix elements). This
%    subfunction is needed when using memory save option in 
%    gp_set.
%
%  See also
%    GPCF_scaled_PAK, GPCF_scaled_UNPAK, GPCF_scaled_LP, GP_G
  
  [n, m] =size(x);
  if nargin==4
    % Use memory save option
    savememory=1;
    if i1==0
      % Return number of covariates
      if isfield(gpcf,'selectedVariables')
        DKff=length(gpcf.selectedVariables);
      else
        DKff=m;
      end
      return
    end
  else
    savememory=0;
  end

  % evaluate the gradient for training covariance
  if nargin == 2 || isempty(x2)
    scale = sparse(1:n,1:n,x(:,gpcf.scaler),n,n);
    DKff = {};
    
    cf = gpcf.cf{1};
    if ~savememory
      DK = cf.fh.ginput(cf, x);
    else
      DK = cf.fh.ginput(cf,x,[],i1);
    end
    
    for j = 1:length(DK)
        DKff{end+1} = scale*DK{j}*scale;
    end

    % Evaluate the gradient of non-symmetric covariance (e.g. K_fu)
  elseif nargin == 3 || nargin == 4
      if size(x,2) ~= size(x2,2)
          error('gpcf_scaled -> _ghyper: The number of columns in x and x2 has to be the same. ')
      end
      scale = sparse(1:n,1:n,x(:,gpcf.scaler),n,n);
      n2 = length(x2);
      scale2 = sparse(1:n2,1:n2,x2(:,gpcf.scaler),n2,n2);
      
      cf = gpcf.cf{1};
      if ~savememory
        DK = cf.fh.ginput(cf, x, x2);
      else
        DK = cf.fh.ginput(cf,x,x2,i1);
      end
      
      for j = 1:length(DK)
          DKff{end+1} = scale*DK{j}*scale2;
      end
  end
  
end


function C = gpcf_scaled_cov(gpcf, x1, x2)
%GP_scaled_COV  Evaluate covariance matrix between two input vectors
%
%  Description         
%    C = GP_scaled_COV(GP, TX, X) takes in covariance function of a
%    Gaussian process GP and two matrixes TX and X that contain
%    input vectors to GP. Returns covariance matrix C. Every
%    element ij of C contains covariance between inputs i in TX
%    and j in X. This is a mandatory subfunction used for example in
%    prediction and energy computations.
%
%
%  See also
%    GPCF_scaled_TRCOV, GPCF_scaled_TRVAR, GP_COV, GP_TRCOV
  
  if isempty(x2)
    x2=x1;
  end
  [n1,m1]=size(x1);
  [n2,m2]=size(x2);
  scale = sparse(1:n1,1:n1,x1(:,gpcf.scaler),n1,n1);
  scale2 = sparse(1:n2,1:n2,x2(:,gpcf.scaler),n2,n2);

  if m1~=m2
    error('the number of columns of X1 and X2 has to be same')
  end

  cf = gpcf.cf{1};
  C = scale*feval(cf.fh.cov, cf, x1, x2)*scale2;
end

function C = gpcf_scaled_trcov(gpcf, x)
%GP_scaled_TRCOV     Evaluate training covariance matrix of inputs
%
%  Description
%    C = GP_scaled_TRCOV(GP, TX) takes in covariance function of a
%    Gaussian process GP and matrix TX that contains training
%    input vectors. Returns covariance matrix C. Every element ij
%    of C contains covariance between inputs i and j in TX. This 
%    is a mandatory subfunction used for example in prediction 
%    and energy computations.
%
%  See also
%    GPCF_scaled_COV, GPCF_scaled_TRVAR, GP_COV, GP_TRCOV
  n = length(x);
  scale = sparse(1:n,1:n,x(:,gpcf.scaler),n,n);
  
  cf = gpcf.cf{1};
  C = scale*feval(cf.fh.trcov, cf, x)*scale;
end

function C = gpcf_scaled_trvar(gpcf, x)
% GP_scaled_TRVAR     Evaluate training variance vector
%
%  Description
%    C = GP_scaled_TRVAR(GPCF, TX) takes in covariance function of
%    a Gaussian process GPCF and matrix TX that contains training
%    inputs. Returns variance vector C. Every element i of C
%    contains variance of input i in TX. This is a mandatory 
%    subfunction used for example in prediction and energy computations.
%
%  See also
%    GPCF_scaled_COV, GP_COV, GP_TRCOV

    cf = gpcf.cf{1};
    C = x(:,gpcf.scaler).*feval(cf.fh.trvar, cf, x).*x(:,gpcf.scaler);
end

function reccf = gpcf_scaled_recappend(reccf, ri, gpcf)
%RECAPPEND  Record append
%
%  Description
%    RECCF = GPCF_scaled_RECAPPEND(RECCF, RI, GPCF) takes a
%    covariance function record structure RECCF, record index RI
%    and covariance function structure GPCF with the current MCMC
%    samples of the parameters. Returns RECCF which contains all
%    the old samples and the current samples from GPCF. This 
%    subfunction is needed when using MCMC sampling (gp_mc).
%
%  See also
%    GP_MC, GP_MC->RECAPPEND
  
% Initialize record
  if nargin == 2
    reccf.type = 'gpcf_scaled';

    cf = ri.cf{1};
    reccf.cf{1} = feval(cf.fh.recappend, [], ri.cf{1});
    reccf.scaler = ri.scaler;
    
    % Set the function handles
    reccf.fh.pak = @gpcf_scaled_pak;
    reccf.fh.unpak = @gpcf_scaled_unpak;
    reccf.fh.e = @gpcf_scaled_lp;
    reccf.fh.lpg = @gpcf_scaled_lpg;
    reccf.fh.cfg = @gpcf_scaled_cfg;
    reccf.fh.cov = @gpcf_scaled_cov;
    reccf.fh.trcov  = @gpcf_scaled_trcov;
    reccf.fh.trvar  = @gpcf_scaled_trvar;
    reccf.fh.recappend = @gpcf_scaled_recappend;
    return
  end
  
  %loop over all of the covariance functions
  cf = gpcf.cf{1};
  reccf.cf{1} = feval(cf.fh.recappend, reccf.cf{1}, ri, cf);
  
end

