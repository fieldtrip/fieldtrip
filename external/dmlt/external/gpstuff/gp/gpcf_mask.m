function gpcf = gpcf_mask(varargin)
%GPCF_MASK  Create a mask covariance function
%
%  Description
%    GPCF = GPCF_MASK('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates a mask covariance function structure in
%    which the named parameters have the specified values. Any
%    unspecified parameters are set to default values.
%
%    GPCF = GPCF_MASK(GPCF,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a covariance function structure with the named
%    parameters altered with the specified values.
%  
%    Mask covariance function returns correlation 1 if input
%    values X_i and X_j are both non zero and 0 otherwise.
%
%    Parameters for mask covariance function
%      selectedVariables - vector defining which inputs are used
%
%  See also
%    GP_SET, GPCF_*, PRIOR_*, MEAN_*

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2008-2010 Jaakko Riihim�ki
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPCF_MASK';
  ip.addOptional('gpcf', [], @isstruct);
  ip.addParamValue('selectedVariables',[], ...
                   @(x) isempty(x) || (isvector(x) && all(x>0)));
  ip.parse(varargin{:});
  gpcf=ip.Results.gpcf;

  if isempty(gpcf)
    % Initialize a covariance function structure
    init=true;
    gpcf.type = 'gpcf_mask';
  else
    % Modify a covariance function structure
    if ~isfield(gpcf,'type') && ~isequal(gpcf.type,'gpcf_mask')
      error('First argument does not seem to be a valid covariance function structure')
    end
    init=false;
  end
  
  if ~ismember('selectedVariables',ip.UsingDefaults)
    if ~isempty(ip.Results.selectedVariables)
      gpcf.selectedVariables = ip.Results.selectedVariables;
    elseif isfield(gpcf,'selectedVariables')
      gpcf=rmfield(gpcf,'selectedVariables');
    end
  end
  if init
    gpcf.fh.pak = @gpcf_mask_pak;
    gpcf.fh.unpak = @gpcf_mask_unpak;
    gpcf.fh.lp = @gpcf_mask_lp;
    gpcf.fh.lpg = @gpcf_mask_lpg;
    gpcf.fh.cfg = @gpcf_mask_cfg;
    gpcf.fh.ginput = @gpcf_mask_ginput;
    gpcf.fh.cov = @gpcf_mask_cov;
    gpcf.fh.trcov  = @gpcf_mask_trcov;
    gpcf.fh.trvar  = @gpcf_mask_trvar;
    gpcf.fh.recappend = @gpcf_mask_recappend;
  end
  
end

function [w,s] = gpcf_mask_pak(gpcf, w)
%GPCF_MASK_PAK  Combine GP covariance function parameters into
%              one vector.
%
%  Description
%    W = GPCF_MASK_PAK(GPCF) takes a covariance function
%    structure GPCF and combines the covariance function
%    parameters and their hyperparameters into a single row
%    vector W. This is a mandatory subfunction used for 
%    example in energy and gradient computations.
%
%       w = []
%
%  See also
%    GPCF_MASK_UNPAK
  
  w = []; s = {};
end

function [gpcf, w] = gpcf_mask_unpak(gpcf, w)
%GPCF_MASK_UNPAK  Sets the covariance function parameters into
%                the structure
%
%  Description
%    [GPCF, W] = GPCF_MASK_UNPAK(GPCF, W) takes a covariance
%    function structure GPCF and a parameter vector W, and
%    returns a covariance function structure identical to the
%    input, except that the covariance parameters have been set
%    to the values in W. Deletes the values set to GPCF from W
%    and returns the modified W. This is a mandatory subfunction 
%    used for example in energy and gradient computations.
%
%    Assignment is inverse of  
%       w = []
%
%  See also
%   GPCF_MASK_PAK
  
end

function lp = gpcf_mask_lp(gpcf)
%GPCF_MASK_LP  Evaluate the energy of prior of covariance function parameters
%
%  Description
%    LP = GPCF_MASK_LP(GPCF) takes a covariance function
%    structure GPCF and returns log(p(th)), where th collects the
%    parameters. This is a mandatory subfunction used for example 
%    in energy computations.
%
%  See also
%    GPCF_MASK_PAK, GPCF_MASK_UNPAK, GPCF_MASK_LPG, GP_E

  lp = 0;
  
end

function lpg = gpcf_mask_lpg(gpcf)
%GPCF_MASK_LPG  Evaluate gradient of the log prior with respect
%               to the parameters.
%
%  Description
%    LPG = GPCF_MASK_LPG(GPCF) takes a covariance function
%    structure GPCF and returns LPG = d log (p(th))/dth, where th
%    is the vector of parameters. This is a mandatory subfunction 
%    used for example in gradient computations.
%
%  See also
%    GPCF_MASK_PAK, GPCF_MASK_UNPAK, GPCF_MASK_LP, GP_G
  
  lpg = [];

end

function DKff = gpcf_mask_cfg(gpcf, x, x2, mask,i1)
%GPCF_MASK_CFG  Evaluate gradient of covariance function
%               with respect to the parameters.
%
%  Description
%    DKff = GPCF_MASK_CFG(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to th (cell array with matrix elements). This is a 
%    mandatory subfunction used in gradient computations.
%
%    DKff = GPCF_MASK_CFG(GPCF, X, X2) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X2) with
%    respect to th (cell array with matrix elements). This subfunction 
%    is needed when using sparse approximations (e.g. FIC).
%
%    DKff = GPCF_MASK_CFG(GPCF, X, [], MASK) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the diagonal of gradients of covariance matrix
%    Kff = k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_MASK_CFG(GPCF, X, X2, [], i) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X2), or k(X,X)
%    if X2 is empty, with respect to ith hyperparameter. This subfunction
%    is needed when using memory save option in gp_set.
%
%  See also
%    GPCF_MASK_PAK, GPCF_MASK_UNPAK, GPCF_MASK_LP, GP_G

  DKff = {};
  
end

function [DKff, lpg]  = gpcf_mask_ginput(gpcf, x, x2, i1)
%GPCF_MASK_GINPUT  Evaluate gradient of covariance function with 
%                  respect to x.
%
%  Description
%    DKff = GPCF_MASK_GINPUT(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to X (cell array with matrix elements). This subfunction
%    is needed when computing gradients with respect to inducing
%    inputs in sparse approximations.
%
%    DKff = GPCF_MASK_GINPUT(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors
%    and returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_MASK_GINPUT(GPCF, X, X2, i) takes a covariance
%    function structure GPCF, a matrix X of input vectors
%    and returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2), or k(X,X) if X2 is empty, with respect to ith covariate
%    in X. This subfunction is needed when using memory saving option
%    in gp_set.
%
%  See also
%   GPCF_MASK_PAK, GPCF_MASK_UNPAK, GPCF_MASK_LP, GP_G
  
  [n, m] =size(x);
  
  if nargin == 2 || isempty(x2)
    ii1 = 0;
    for i=1:m
      for j = 1:n
        ii1 = ii1 + 1;
        DKff{ii1} = zeros(n);
        lpg(ii1) = 0;
      end
    end
  elseif nargin == 3 || nargin == 4
    ii1 = 0;
    for i=1:m
      for j = 1:n
        ii1 = ii1 + 1;
        DKff{ii1} = zeros(n, size(x2,1));
        lpg(ii1) = 0; 
      end
    end
  end
  if nargin==4
    DKff=DKff{1};
  end
end

function C = gpcf_mask_cov(gpcf, x1, x2, varargin)
%GP_MASK_COV  Evaluate covariance matrix between two input vectors
%
%  Description         
%    C = GP_MASK_COV(GP, TX, X) takes in covariance function of a
%    Gaussian process GP and two matrixes TX and X that contain
%    input vectors to GP. Returns covariance matrix C. Every
%    element ij of C contains covariance between inputs i in TX
%    and j in X. This is a mandatory subfunction used for example in
%    prediction and energy computations.
%
%  See also
%    GPCF_MASK_TRCOV, GPCF_MASK_TRVAR, GP_COV, GP_TRCOV
  
  if isempty(x2)
    x2=x1;
  end
  [n1,m1]=size(x1);
  [n2,m2]=size(x2);
  
  if m1~=m2
    error('the number of columns of X1 and X2 has to be same')
  end

  C=repmat(true,n1,n2);
  if isfield(gpcf, 'selectedVariables')
    for j = 1:length(gpcf.selectedVariables)
      jj=gpcf.selectedVariables(j);
      C = C & bsxfun(@and,x1(:,jj)~=0,x2(:,jj)'~=0);
    end
  else
    for j = 1:m1
      C = C & bsxfun(@and,x1(:,j)~=0,x2(:,j)'~=0);
    end
  end
  C=double(C);
  
end

function C = gpcf_mask_trcov(gpcf, x)
%GP_MASK_TRCOV  Evaluate training covariance matrix of inputs
%
%  Description
%    C = GP_MASK_TRCOV(GP, TX) takes in covariance function of a
%    Gaussian process GP and matrix TX that contains training
%    input vectors. Returns covariance matrix C. Every element ij
%    of C contains covariance between inputs i and j in TX. This 
%    is a mandatory subfunction used for example in prediction 
%    and energy computations.
%
%  See also
%    GPCF_MASK_COV, GPCF_MASK_TRVAR, GP_COV, GP_TRCOV

  [n,m]=size(x);

  C=repmat(true,n,n);
  if isfield(gpcf, 'selectedVariables')
    for j = 1:length(gpcf.selectedVariables)
      jj=gpcf.selectedVariables(j);
      C = C & bsxfun(@and,x(:,jj)~=0,x(:,jj)'~=0);
    end
  else
    for j = 1:m
      C = C & bsxfun(@and,x(:,j)~=0,x(:,j)'~=0);
    end
  end
  C=double(C);
  
end

function C = gpcf_mask_trvar(gpcf, x)
%GP_MASK_TRVAR  Evaluate training variance vector
%
%  Description
%    C = GP_MASK_TRVAR(GPCF, TX) takes in covariance function of a
%    Gaussian process GPCF and matrix TX that contains training
%    inputs. Returns variance vector C. Every element i of C
%    contains variance of input i in TX. This is a mandatory 
%    subfunction used for example in prediction and energy 
%    computations.
%
%  See also
%    GPCF_MASK_COV, GP_COV, GP_TRCOV

  [n,m]=size(x);
  C=true(n,1);
  if isfield(gpcf, 'selectedVariables')
    for j = 1:length(gpcf.selectedVariables)
      jj=gpcf.selectedVariables(j);
      C = C & x(:,jj)~=0;
    end
  else
    for j = 1:m
      C = C & x(:,j)~=0;
    end
  end
  C=double(C);
  
end

function reccf = gpcf_mask_recappend(reccf, ri, gpcf)
%RECAPPEND  Record append
%
%  Description
%    RECCF = GPCF_MASK_RECAPPEND(RECCF, RI, GPCF) takes a
%    covariance function record structure RECCF, record index RI
%    and covariance function structure GPCF with the current MCMC
%    samples of the parameters. Returns RECCF which contains
%    all the old samples and the current samples from GPCF.
%    This subfunction is needed when using MCMC sampling (gp_mc).
%
%  See also
%    GP_MC and GP_MC -> RECAPPEND

  if nargin == 2
    % Initialize the record
    reccf.type = 'gpcf_mask';

    % Initialize parameters
    reccf.coeffSigma2= [];

    % Set the function handles
    reccf.fh.pak = @gpcf_mask_pak;
    reccf.fh.unpak = @gpcf_mask_unpak;
    reccf.fh.lp = @gpcf_mask_lp;
    reccf.fh.lpg = @gpcf_mask_lpg;
    reccf.fh.cfg = @gpcf_mask_cfg;
    reccf.fh.cov = @gpcf_mask_cov;
    reccf.fh.trcov  = @gpcf_mask_trcov;
    reccf.fh.trvar  = @gpcf_mask_trvar;
    reccf.fh.recappend = @gpcf_mask_recappend;
  else
    % Append to the record
    if isfield(gpcf, 'selectedVariables')
      reccf.selectedVariables = gpcf.selectedVariables;
    end
  end
end
