function gpcf = gpcf_noise(varargin)
%GPCF_NOISE  Create a independent noise covariance function
%
%  Description
%    GPCF = GPCF_NOISE('PARAM1',VALUE1,'PARAM2,VALUE2,...) creates
%    independent noise covariance function structure in which the
%    named parameters have the specified values. Any unspecified
%    parameters are set to default values.
%
%    GPCF = GPCF_NOISE(GPCF,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a covariance function structure with the named
%    parameters altered with the specified values.
%
%    Parameters for independent noise covariance function [default]
%      noiseSigma2       - variance of the independent noise [0.1]
%      noiseSigma2_prior - prior for noiseSigma2 [prior_logunif]
%
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc.
%
%  See also
%    GP_SET, GPCF_*, PRIOR_*

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPCF_NOISE';
  ip.addOptional('gpcf', [], @isstruct);
  ip.addParamValue('noiseSigma2',0.1, @(x) isscalar(x) && x>0);
  ip.addParamValue('noiseSigma2_prior',prior_logunif, @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  gpcf=ip.Results.gpcf;
  
  if isempty(gpcf)
    init=true;
    gpcf.type = 'gpcf_noise';
  else
    if ~isfield(gpcf,'type') && ~isequal(gpcf.type,'gpcf_noise')
      error('First argument does not seem to be a valid covariance function structure')
    end
    init=false;
  end
  
  % Initialize parameter
  if init || ~ismember('noiseSigma2',ip.UsingDefaults)
    gpcf.noiseSigma2=ip.Results.noiseSigma2;
  end
  
  % Initialize prior structure
  if init
    gpcf.p=[];
  end
  if init || ~ismember('noiseSigma2_prior',ip.UsingDefaults)
    gpcf.p.noiseSigma2=ip.Results.noiseSigma2_prior;
  end
  
  if init
    % Set the function handles to the subfunctions
    gpcf.fh.pak = @gpcf_noise_pak;
    gpcf.fh.unpak = @gpcf_noise_unpak;
    gpcf.fh.lp = @gpcf_noise_lp;
    gpcf.fh.lpg = @gpcf_noise_lpg;
    gpcf.fh.cfg = @gpcf_noise_cfg;
    gpcf.fh.ginput = @gpcf_noise_ginput;
    gpcf.fh.cov = @gpcf_noise_cov;
    gpcf.fh.trcov  = @gpcf_noise_trcov;
    gpcf.fh.trvar  = @gpcf_noise_trvar;
    gpcf.fh.recappend = @gpcf_noise_recappend;
  end        

end

function [w, s] = gpcf_noise_pak(gpcf)
%GPCF_NOISE_PAK  Combine GP covariance function parameters into
%                one vector.
%
%  Description
%    W = GPCF_NOISE_PAK(GPCF) takes a covariance function data
%    structure GPCF and combines the covariance function
%    parameters and their hyperparameters into a single row
%    vector W. This is a mandatory subfunction used for example 
%    in energy and gradient computations.
%
%       w = [ log(gpcf.noiseSigma2)
%             (hyperparameters of gpcf.magnSigma2)]'
%     
%
%  See also
%    GPCF_NOISE_UNPAK


  w = []; s = {};
  if ~isempty(gpcf.p.noiseSigma2)
    w(1) = log(gpcf.noiseSigma2);
    s = [s 'log(noise.noiseSigma2)'];
    % Hyperparameters of noiseSigma2
    [wh sh] = gpcf.p.noiseSigma2.fh.pak(gpcf.p.noiseSigma2);
    w = [w wh];
    s = [s sh];
  end    

end

function [gpcf, w] = gpcf_noise_unpak(gpcf, w)
%GPCF_NOISE_UNPAK  Sets the covariance function parameters
%                     into the structure
%
%  Description
%    [GPCF, W] = GPCF_NOISE_UNPAK(GPCF, W) takes a covariance
%    function data structure GPCF and a hyper-parameter vector W,
%    and returns a covariance function data structure identical
%    to the input, except that the covariance hyper-parameters
%    have been set to the values in W. Deletes the values set to
%    GPCF from W and returns the modified W. This is a mandatory 
%    subfunction used for example in energy and gradient computations.
%
%    Assignment is inverse of  
%       w = [ log(gpcf.noiseSigma2)
%             (hyperparameters of gpcf.magnSigma2)]'
%
%  See also
%    GPCF_NOISE_PAK

  
  if ~isempty(gpcf.p.noiseSigma2)
    gpcf.noiseSigma2 = exp(w(1));
    w = w(2:end);
    
    % Hyperparameters of lengthScale
    [p, w] = gpcf.p.noiseSigma2.fh.unpak(gpcf.p.noiseSigma2, w);
    gpcf.p.noiseSigma2 = p;
  end
end


function lp = gpcf_noise_lp(gpcf)
%GPCF_NOISE_LP  Evaluate the log prior of covariance function parameters
%
%  Description
%    LP = GPCF_NOISE_LP(GPCF) takes a covariance function
%    structure GPCF and returns log(p(th)), where th collects the
%    parameters. This is a mandatory subfunction used for example 
%    in energy computations.
%
%  See also
%   GPCF_NOISE_PAK, GPCF_NOISE_UNPAK, GPCF_NOISE_G, GP_E

% Evaluate the prior contribution to the error. The parameters that
% are sampled are from space W = log(w) where w is all the
% "real" samples. On the other hand errors are evaluated in the
% W-space so we need take into account also the Jacobian of
% transformation W -> w = exp(W). See Gelman et.al., 2004,
% Bayesian data Analysis, second edition, p24.
  
  lp = 0;
  gpp=gpcf.p;

  if ~isempty(gpcf.p.noiseSigma2)
    % Evaluate the prior contribution to the error.
    lp = gpp.noiseSigma2.fh.lp(gpcf.noiseSigma2, gpp.noiseSigma2) +log(gpcf.noiseSigma2);
  end
end

function lpg = gpcf_noise_lpg(gpcf)
%GPCF_NOISE_LPG  Evaluate gradient of the log prior with respect
%                to the parameters.
%
%  Description
%    LPG = GPCF_NOISE_LPG(GPCF) takes a covariance function
%    structure GPCF and returns LPG = d log (p(th))/dth, where th
%    is the vector of parameters. This is a mandatory subfunction 
%    used for example in gradient computations.
%
%  See also
%    GPCF_NOISE_PAK, GPCF_NOISE_UNPAK, GPCF_NOISE_LP, GP_G

  lpg = [];
  gpp=gpcf.p;
  
  if ~isempty(gpcf.p.noiseSigma2)            
    lpgs = gpp.noiseSigma2.fh.lpg(gpcf.noiseSigma2, gpp.noiseSigma2);
    lpg = [lpg lpgs(1).*gpcf.noiseSigma2+1 lpgs(2:end)];
  end
end

function DKff = gpcf_noise_cfg(gpcf, x, x2, mask, i1)
%GPCF_NOISE_CFG  Evaluate gradient of covariance function
%                with respect to the parameters
%
%  Description
%    DKff = GPCF_NOISE_CFG(GPCF, X) takes a covariance function
%    data structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to th (cell array with matrix elements). This is a 
%    mandatory subfunction used in gradient computations.
%
%    DKff = GPCF_NOISE_CFG(GPCF, X, X2) takes a covariance
%    function data structure GPCF, a matrix X of input vectors
%    and returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_NOISE_CFG(GPCF, X, [], MASK) takes a covariance
%    function data structure GPCF, a matrix X of input vectors
%    and returns DKff, the diagonal of gradients of covariance
%    matrix Kff = k(X,X2) with respect to th (cell array with
%    matrix elements). This subfunction is needed when using 
%    sparse approximations (e.g. FIC).
%
%  See also
%    GPCF_NOISE_PAK, GPCF_NOISE_UNPAK, GPCF_NOISE_E, GP_G

  DKff = {};

  if ~isempty(gpcf.p.noiseSigma2)
    gpp=gpcf.p;
    DKff{1}=gpcf.noiseSigma2;
  end
  if nargin==4
    % Use memory save option
    if i1==0
      % Return number of hyperparameters
      DKff=1;
      return
    end
    DKff=DKff{1};
  end
  
end

function DKff = gpcf_noise_ginput(gpcf, x, t, i1)
%GPCF_NOISE_GINPUT  Evaluate gradient of covariance function with 
%                   respect to x
%
%  Description
%    DKff = GPCF_NOISE_GINPUT(GPCF, X) takes a covariance
%    function data structure GPCF, a matrix X of input vectors
%    and returns DKff, the gradients of covariance matrix Kff =
%    k(X,X) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_NOISE_GINPUT(GPCF, X, X2) takes a covariance
%    function data structure GPCF, a matrix X of input vectors
%    and returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%  See also
%    GPCF_NOISE_PAK, GPCF_NOISE_UNPAK, GPCF_NOISE_E, GP_G

end

function C = gpcf_noise_cov(gpcf, x1, x2)
% GP_NOISE_COV     Evaluate covariance matrix between two input vectors
%
%  Description
%    C = GP_NOISE_COV(GP, TX, X) takes in covariance function of
%    a Gaussian process GP and two matrixes TX and X that contain
%    input vectors to GP. Returns covariance matrix C. Every
%    element ij of C contains covariance between inputs i in TX
%    and j in X. This is a mandatory subfunction used for example in
%    prediction and energy computations.
%
%  See also
%    GPCF_NOISE_TRCOV, GPCF_NOISE_TRVAR, GP_COV, GP_TRCOV


  if isempty(x2)
    x2=x1;
  end
  [n1,m1]=size(x1);
  [n2,m2]=size(x2);

  if m1~=m2
    error('the number of columns of X1 and X2 has to be same')
  end

  C = sparse([],[],[],n1,n2,0);
end

function C = gpcf_noise_trcov(gpcf, x)
%GP_NOISE_TRCOV  Evaluate training covariance matrix of inputs
%
%  Description
%    C = GP_NOISE_TRCOV(GP, TX) takes in covariance function of a
%    Gaussian process GP and matrix TX that contains training
%    input vectors. Returns covariance matrix C. Every element ij
%    of C contains covariance between inputs i and j in TX. This is 
%    a mandatory subfunction used for example in prediction and 
%    energy computations.
%   
%
%  See also
%    GPCF_NOISE_COV, GPCF_NOISE_TRVAR, GP_COV, GP_TRCOV

  [n, m] =size(x);
  n1=n+1;

  C = sparse([],[],[],n,n,0);
  C(1:n1:end)=C(1:n1:end)+gpcf.noiseSigma2;

end

function C = gpcf_noise_trvar(gpcf, x)
% GP_NOISE_TRVAR     Evaluate training variance vector
%
%    Description
%    C = GP_NOISE_TRVAR(GPCF, TX) takes in covariance function 
%    of a Gaussian process GPCF and matrix TX that contains
%    training inputs. Returns variance vector C. Every
%    element i of C contains variance of input i in TX. This is 
%    a mandatory subfunction used for example in prediction and 
%    energy computations.
%
%
%    See also
%    GPCF_NOISE_COV, GP_COV, GP_TRCOV



  [n, m] =size(x);
  C=ones(n,1)*gpcf.noiseSigma2;

end

function reccf = gpcf_noise_recappend(reccf, ri, gpcf)
%RECAPPEND Record append
%
%  Description
%    RECCF = GPCF_NOISE_RECAPPEND(RECCF, RI, GPCF) takes a
%    covariance function record structure RECCF, record index RI
%    and covariance function structure GPCF with the current MCMC
%    samples of the hyperparameters. Returns RECCF which contains
%    all the old samples and the current samples from GPCF.
%    This subfunction is needed when using MCMC sampling (gp_mc).
%
%  See also
%    GP_MC and GP_MC -> RECAPPEND

  if nargin == 2
    % Initialize the record
    reccf.type = 'gpcf_noise';
    
    % Initialize parameters
    reccf.noiseSigma2 = []; 
    
    % Set the function handles
    reccf.fh.pak = @gpcf_noise_pak;
    reccf.fh.unpak = @gpcf_noise_unpak;
    reccf.fh.e = @gpcf_noise_lp;
    reccf.fh.lpg = @gpcf_noise_lpg;
    reccf.fh.cfg = @gpcf_noise_cfg;
    reccf.fh.cov = @gpcf_noise_cov;
    reccf.fh.trcov  = @gpcf_noise_trcov;
    reccf.fh.trvar  = @gpcf_noise_trvar;
    %  gpcf.fh.sampling = @hmc2;
    reccf.sampling_opt = hmc2_opt;
    reccf.fh.recappend = @gpcf_noise_recappend;  
    reccf.p=[];
    reccf.p.noiseSigma2=[];
    if ~isempty(ri.p.noiseSigma2)
      reccf.p.noiseSigma2 = ri.p.noiseSigma2;
    end
  else
    % Append to the record
    gpp = gpcf.p;

    % record noiseSigma2
    reccf.noiseSigma2(ri,:)=gpcf.noiseSigma2;
    if isfield(gpp,'noiseSigma2') && ~isempty(gpp.noiseSigma2)
      reccf.p.noiseSigma2 = gpp.noiseSigma2.fh.recappend(reccf.p.noiseSigma2, ri, gpcf.p.noiseSigma2);
    end
  end
end
