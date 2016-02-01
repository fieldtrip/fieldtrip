function gpcf = gpcf_linear(varargin)
%GPCF_LINEAR  Create a linear (dot product) covariance function
%
%  Description
%    GPCF = GPCF_LINEAR('PARAM1',VALUE1,'PARAM2,VALUE2,...) creates
%    a linear (dot product) covariance function structure in which
%    the named parameters have the specified values. Any
%    unspecified parameters are set to default values.
%
%    GPCF = GPCF_LINEAR(GPCF,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a covariance function structure with the named
%    parameters altered with the specified values.
%  
%    Parameters for linear (dot product) covariance function
%      coeffSigma2       - prior variance for regressor coefficients [10]
%                          This can be either scalar corresponding
%                          to a common prior variance or vector
%                          defining own prior variance for each
%                          coefficient.
%      coeffSigma2_prior - prior structure for coeffSigma2 [prior_logunif]
%      selectedVariables - vector defining which inputs are used [all]
%
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc.
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
  ip.FunctionName = 'GPCF_LINEAR';
  ip.addOptional('gpcf', [], @isstruct);
  ip.addParamValue('coeffSigma2',10, @(x) isvector(x) && all(x>0));
  ip.addParamValue('coeffSigma2_prior',prior_logunif, @(x) isstruct(x) || isempty(x));
  ip.addParamValue('selectedVariables',[], @(x) isvector(x) && all(x>0));
  ip.parse(varargin{:});
  gpcf=ip.Results.gpcf;

  if isempty(gpcf)
    init=true;
    gpcf.type = 'gpcf_linear';
  else
    if ~isfield(gpcf,'type') && ~isequal(gpcf.type,'gpcf_linear')
      error('First argument does not seem to be a valid covariance function structure')
    end
    init=false;
  end
  
  % Initialize parameter
  if init || ~ismember('coeffSigma2',ip.UsingDefaults)
    gpcf.coeffSigma2=ip.Results.coeffSigma2;
  end

  % Initialize prior structure
  if init
    gpcf.p=[];
  end
  if init || ~ismember('coeffSigma2_prior',ip.UsingDefaults)
    gpcf.p.coeffSigma2=ip.Results.coeffSigma2_prior;
  end
  if ~ismember('selectedVariables',ip.UsingDefaults)
    selectedVariables=ip.Results.selectedVariables;
    if ~isempty(selectedVariables)
      gpcf.selectedVariables = selectedVariables;
    end
  end
  
  if init
    % Set the function handles to the subfunctions
    gpcf.fh.pak = @gpcf_linear_pak;
    gpcf.fh.unpak = @gpcf_linear_unpak;
    gpcf.fh.lp = @gpcf_linear_lp;
    gpcf.fh.lpg = @gpcf_linear_lpg;
    gpcf.fh.cfg = @gpcf_linear_cfg;
    gpcf.fh.ginput = @gpcf_linear_ginput;
    gpcf.fh.cov = @gpcf_linear_cov;
    gpcf.fh.trcov  = @gpcf_linear_trcov;
    gpcf.fh.trvar  = @gpcf_linear_trvar;
    gpcf.fh.recappend = @gpcf_linear_recappend;
  end        

end

function [w, s] = gpcf_linear_pak(gpcf, w)
%GPCF_LINEAR_PAK  Combine GP covariance function parameters into one vector
%
%  Description
%    W = GPCF_LINEAR_PAK(GPCF) takes a covariance function
%    structure GPCF and combines the covariance function
%    parameters and their hyperparameters into a single row
%    vector W. This is a mandatory subfunction used for 
%    example in energy and gradient computations.
%
%       w = [ log(gpcf.coeffSigma2)
%             (hyperparameters of gpcf.coeffSigma2)]'
%
%  See also
%    GPCF_LINEAR_UNPAK
  
  w = []; s = {};
  if ~isempty(gpcf.p.coeffSigma2)
    w = log(gpcf.coeffSigma2);
    if numel(gpcf.coeffSigma2)>1
      s = [s; sprintf('log(linear.coeffSigma2 x %d)',numel(gpcf.coeffSigma2))];
    else
      s = [s; 'log(linear.coeffSigma2)'];
    end
    % Hyperparameters of coeffSigma2
    [wh sh] = gpcf.p.coeffSigma2.fh.pak(gpcf.p.coeffSigma2);
    w = [w wh];
    s = [s; sh];
  end
end

function [gpcf, w] = gpcf_linear_unpak(gpcf, w)
%GPCF_LINEAR_UNPAK  Sets the covariance function parameters 
%                   into the structure
%
%  Description
%    [GPCF, W] = GPCF_LINEAR_UNPAK(GPCF, W) takes a covariance
%    function structure GPCF and a hyper-parameter vector W, and
%    returns a covariance function structure identical to the
%    input, except that the covariance hyper-parameters have been
%    set to the values in W. Deletes the values set to GPCF from
%    W and returns the modified W. This is a mandatory subfunction 
%    used for example in energy and gradient computations.
%
%    Assignment is inverse of  
%       w = [ log(gpcf.coeffSigma2)
%             (hyperparameters of gpcf.coeffSigma2)]'
%
%  See also
%   GPCF_LINEAR_PAK
  
  gpp=gpcf.p;

  if ~isempty(gpp.coeffSigma2)
    i2=length(gpcf.coeffSigma2);
    i1=1;
    gpcf.coeffSigma2 = exp(w(i1:i2));
    w = w(i2+1:end);
    
    % Hyperparameters of coeffSigma2
    [p, w] = gpcf.p.coeffSigma2.fh.unpak(gpcf.p.coeffSigma2, w);
    gpcf.p.coeffSigma2 = p;
  end
end

function lp = gpcf_linear_lp(gpcf)
%GPCF_LINEAR_LP  Evaluate the log prior of covariance function parameters
%
%  Description
%    LP = GPCF_LINEAR_LP(GPCF) takes a covariance function
%    structure GPCF and returns log(p(th)), where th collects the
%    parameters. This is a mandatory subfunction used for example 
%    in energy computations.
%
%  See also
%   GPCF_LINEAR_PAK, GPCF_LINEAR_UNPAK, GPCF_LINEAR_LPG, GP_E

% Evaluate the prior contribution to the error. The parameters that
% are sampled are from space W = log(w) where w is all the "real" samples.
% On the other hand errors are evaluated in the W-space so we need take
% into account also the  Jacobian of transformation W -> w = exp(W).
% See Gelman et.al., 2004, Bayesian data Analysis, second edition, p24.
  lp = 0;
  gpp=gpcf.p;

  if ~isempty(gpp.coeffSigma2)
    lp = gpp.coeffSigma2.fh.lp(gpcf.coeffSigma2, gpp.coeffSigma2) + sum(log(gpcf.coeffSigma2));
  end
end

function lpg = gpcf_linear_lpg(gpcf)
%GPCF_LINEAR_LPG  Evaluate gradient of the log prior with respect
%                 to the parameters.
%
%  Description
%    LPG = GPCF_LINEAR_LPG(GPCF) takes a covariance function
%    structure GPCF and returns LPG = d log (p(th))/dth, where th
%    is the vector of parameters. This is a mandatory subfunction 
%    used for example in gradient computations.
%
%  See also
%    GPCF_LINEAR_PAK, GPCF_LINEAR_UNPAK, GPCF_LINEAR_LP, GP_G

  lpg = [];
  gpp=gpcf.p;
  
  if ~isempty(gpcf.p.coeffSigma2)            
    lll=length(gpcf.coeffSigma2);
    lpgs = gpp.coeffSigma2.fh.lpg(gpcf.coeffSigma2, gpp.coeffSigma2);
    lpg = [lpg lpgs(1:lll).*gpcf.coeffSigma2+1 lpgs(lll+1:end)];
  end
end

function DKff = gpcf_linear_cfg(gpcf, x, x2, mask, i1)
%GPCF_LINEAR_CFG  Evaluate gradient of covariance function
%                 with respect to the parameters
%
%  Description
%    DKff = GPCF_LINEAR_CFG(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to th (cell array with matrix elements). This is a 
%    mandatory subfunction used in gradient computations.
%
%    DKff = GPCF_LINEAR_CFG(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_LINEAR_CFG(GPCF, X, [], MASK) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the diagonal of gradients of covariance matrix
%    Kff = k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_LINEAR_CFG(GPCF,X,X2,MASK,i) takes a covariance 
%    function structure GPCF, a matrix X of input vectors and 
%    returns DKff, the gradient of covariance matrix Kff = 
%    k(X,X2), or k(X,X) if X2 is empty, with respect to ith 
%    hyperparameter. This subfunction is needed when using
%    memory save option in gp_set.
%
%  See also
%   GPCF_LINEAR_PAK, GPCF_LINEAR_UNPAK, GPCF_LINEAR_LP, GP_G

  [n, m] =size(x);

  DKff = {};
  
  if nargin==5
    % Use memory save option
    savememory=1;
    if i1==0
      % Return number of hyperparameters
      DKff=0;
      if ~isempty(gpcf.p.coeffSigma2)
        DKff=length(gpcf.coeffSigma2);
      end
      return
    end
  else
    savememory=0;
  end
  
  % Evaluate: DKff{1} = d Kff / d coeffSigma2
  % NOTE! Here we have already taken into account that the parameters are transformed
  % through log() and thus dK/dlog(p) = p * dK/dp

  
  % evaluate the gradient for training covariance
  if nargin == 2 || (isempty(x2) && isempty(mask))
    
    if isfield(gpcf, 'selectedVariables')
      if ~isempty(gpcf.p.coeffSigma2)
        if length(gpcf.coeffSigma2) == 1
          DKff{1}=gpcf.coeffSigma2*x(:,gpcf.selectedVariables)*(x(:,gpcf.selectedVariables)');
        else
          if ~savememory
            i1=1:length(gpcf.coeffSigma2);
          end
          for ii1=i1
            DD = gpcf.coeffSigma2(ii1)*x(:,gpcf.selectedVariables(ii1))*(x(:,gpcf.selectedVariables(ii1))');
            DD(abs(DD)<=eps) = 0;
            DKff{ii1}= (DD+DD')./2;
          end
        end
      end
    else
      if ~isempty(gpcf.p.coeffSigma2)
        if length(gpcf.coeffSigma2) == 1
          DKff{1}=gpcf.coeffSigma2*x*(x');
        else
          if isa(gpcf.coeffSigma2,'single')
            epsi=eps('single');
          else
            epsi=eps;
          end
          if ~savememory
            i1=1:length(gpcf.coeffSigma2);
          end
          DKff=cell(1,length(i1));
          for ii1=i1
            DD = gpcf.coeffSigma2(ii1)*x(:,ii1)*(x(:,ii1)');
            DD(abs(DD)<=epsi) = 0;
            DKff{ii1}= (DD+DD')./2;
          end
        end
      end
    end
    
    
    % Evaluate the gradient of non-symmetric covariance (e.g. K_fu)
  elseif nargin == 3 || isempty(mask)
    if size(x,2) ~= size(x2,2)
      error('gpcf_linear -> _ghyper: The number of columns in x and x2 has to be the same. ')
    end

    if isfield(gpcf, 'selectedVariables')
      if ~isempty(gpcf.p.coeffSigma2)
        if length(gpcf.coeffSigma2) == 1
          DKff{1}=gpcf.coeffSigma2*x(:,gpcf.selectedVariables)*(x2(:,gpcf.selectedVariables)');
        else
          if ~savememory
            i1=1:length(gpcf.coeffSigma2);
          end
          for ii1=i1
            DKff{ii1}=gpcf.coeffSigma2(ii1)*x(:,gpcf.selectedVariables(ii1))*(x2(:,gpcf.selectedVariables(ii1))');
          end
        end
      end
    else
      if ~isempty(gpcf.p.coeffSigma2)
        if length(gpcf.coeffSigma2) == 1
          DKff{1}=gpcf.coeffSigma2*x*(x2');
        else
          if ~savememory
            i1=1:m;
          end            
          for ii1=i1
            DKff{ii1}=gpcf.coeffSigma2(ii1)*x(:,ii1)*(x2(:,ii1)');
          end
        end
      end
    end
    % Evaluate: DKff{1}    = d mask(Kff,I) / d coeffSigma2
    %           DKff{2...} = d mask(Kff,I) / d coeffSigma2
  elseif nargin == 4 || nargin == 5
    
    if isfield(gpcf, 'selectedVariables')
      if ~isempty(gpcf.p.coeffSigma2)
        if length(gpcf.coeffSigma2) == 1
          DKff{1}=gpcf.coeffSigma2*sum(x(:,gpcf.selectedVariables).^2,2); % d mask(Kff,I) / d coeffSigma2
        else
          if ~savememory
            i1=1:length(gpcf.coeffSigma2);
          end
          for ii1=i1
            DKff{ii1}=gpcf.coeffSigma2(ii1)*(x(:,gpcf.selectedVariables(ii1)).^2); % d mask(Kff,I) / d coeffSigma2
          end
        end
      end
    else
      if ~isempty(gpcf.p.coeffSigma2)
        if length(gpcf.coeffSigma2) == 1
          DKff{1}=gpcf.coeffSigma2*sum(x.^2,2); % d mask(Kff,I) / d coeffSigma2
        else
          if ~savememory
            i1=1:m;
          end
          for ii1=i1
            DKff{ii1}=gpcf.coeffSigma2(ii1)*(x(:,ii1).^2); % d mask(Kff,I) / d coeffSigma2
          end
        end
      end
    end
  end
  if savememory
    DKff=DKff{i1};
  end
end


function DKff = gpcf_linear_ginput(gpcf, x, x2, i1)
%GPCF_LINEAR_GINPUT  Evaluate gradient of covariance function with 
%                    respect to x.
%
%  Description
%    DKff = GPCF_LINEAR_GINPUT(GPCF, X) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_LINEAR_GINPUT(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_LINEAR_GINPUT(GPCF, X, X2, i) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to ith covariate in X (matrix).
%    This subfunction is needed when using memory save option
%    in gp_set.
%
%  See also
%   GPCF_LINEAR_PAK, GPCF_LINEAR_UNPAK, GPCF_LINEAR_LP, GP_G        
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
  
  if nargin == 2 || isempty(x2)
    
    %K = feval(gpcf.fh.trcov, gpcf, x);
    
    if length(gpcf.coeffSigma2) == 1
      % In the case of an isotropic LINEAR
      s = repmat(gpcf.coeffSigma2, 1, m);
    else
      s = gpcf.coeffSigma2;
    end
    
    ii1 = 0;
    if isfield(gpcf, 'selectedVariables')
      if ~savememory
        i1=1:length(gpcf.selectedVariables);
      end
      for i=i1
        for j = 1:n
          
          DK = zeros(n);
          DK(j,:)=s(i)*x(:,gpcf.selectedVariables(i))';
          
          DK = DK + DK';
          
          ii1 = ii1 + 1;
          DKff{ii1} = DK;
        end
      end
    else
      if ~savememory
        i1=1:m;
      end
      for i=i1
        for j = 1:n
          
          DK = zeros(n);
          DK(j,:)=s(i)*x(:,i)';
          
          DK = DK + DK';
          
          ii1 = ii1 + 1;
          DKff{ii1} = DK;
        end
      end
    end
    
    
    
  elseif nargin == 3 || nargin == 4
    %K = feval(gpcf.fh.cov, gpcf, x, x2);
    
    if length(gpcf.coeffSigma2) == 1
      % In the case of an isotropic LINEAR
      s = repmat(gpcf.coeffSigma2, 1, m);
    else
      s = gpcf.coeffSigma2;
    end
    
    ii1 = 0;
    if isfield(gpcf, 'selectedVariables')
      if ~savememory
        i1=1:length(gpcf.selectedVariables);
      end
      for i=i1
        for j = 1:n
          
          DK = zeros(n, size(x2,1));
          DK(j,:)=s(i)*x2(:,gpcf.selectedVariables(i))';
          
          ii1 = ii1 + 1;
          DKff{ii1} = DK;
        end
      end
    else
      if ~savememory
        i1=1:m;
      end
      for i=i1
        for j = 1:n
          
          DK = zeros(n, size(x2,1));
          DK(j,:)=s(i)*x2(:,i)';
          
          ii1 = ii1 + 1;
          DKff{ii1} = DK;
        end
      end
    end
    
  end
end


function C = gpcf_linear_cov(gpcf, x1, x2, varargin)
%GP_LINEAR_COV  Evaluate covariance matrix between two input vectors
%
%  Description         
%    C = GP_LINEAR_COV(GP, TX, X) takes in covariance function of
%    a Gaussian process GP and two matrixes TX and X that contain
%    input vectors to GP. Returns covariance matrix C. Every
%    element ij of C contains covariance between inputs i in TX
%    and j in X. This is a mandatory subfunction used for example in
%    prediction and energy computations.
%
%  See also
%    GPCF_LINEAR_TRCOV, GPCF_LINEAR_TRVAR, GP_COV, GP_TRCOV
  
  if isempty(x2)
    x2=x1;
  end
  [n1,m1]=size(x1);
  [n2,m2]=size(x2);

  if m1~=m2
    error('the number of columns of X1 and X2 has to be same')
  end
  
  if isfield(gpcf, 'selectedVariables')
    C = x1(:,gpcf.selectedVariables)*diag(gpcf.coeffSigma2)*(x2(:,gpcf.selectedVariables)');
  else
    C = x1*diag(gpcf.coeffSigma2)*(x2');
  end
  C(abs(C)<=eps) = 0;
end

function C = gpcf_linear_trcov(gpcf, x)
%GP_LINEAR_TRCOV  Evaluate training covariance matrix of inputs
%
%  Description
%    C = GP_LINEAR_TRCOV(GP, TX) takes in covariance function of
%    a Gaussian process GP and matrix TX that contains training
%    input vectors. Returns covariance matrix C. Every element ij
%    of C contains covariance between inputs i and j in TX. This 
%    is a mandatory subfunction used for example in prediction and 
%    energy computations.
%
%  See also
%    GPCF_LINEAR_COV, GPCF_LINEAR_TRVAR, GP_COV, GP_TRCOV

  if isfield(gpcf, 'selectedVariables')
    C = x(:,gpcf.selectedVariables)*diag(gpcf.coeffSigma2)*(x(:,gpcf.selectedVariables)');
  else
    C = x*diag(gpcf.coeffSigma2)*(x');
  end
  C(abs(C)<=eps) = 0;
  C = (C+C')./2;

end


function C = gpcf_linear_trvar(gpcf, x)
%GP_LINEAR_TRVAR  Evaluate training variance vector
%
%  Description
%    C = GP_LINEAR_TRVAR(GPCF, TX) takes in covariance function
%    of a Gaussian process GPCF and matrix TX that contains
%    training inputs. Returns variance vector C. Every element i
%    of C contains variance of input i in TX. This is a mandatory 
%    subfunction used for example in prediction and energy computations.
%
%
%  See also
%    GPCF_LINEAR_COV, GP_COV, GP_TRCOV

  if length(gpcf.coeffSigma2) == 1
    if isfield(gpcf, 'selectedVariables')
      C=gpcf.coeffSigma2.*sum(x(:,gpcf.selectedVariables).^2,2);
    else
      C=gpcf.coeffSigma2.*sum(x.^2,2);
    end
  else
    if isfield(gpcf, 'selectedVariables')
      C=sum(repmat(gpcf.coeffSigma2, size(x,1), 1).*x(:,gpcf.selectedVariables).^2,2);
    else
      C=sum(repmat(gpcf.coeffSigma2, size(x,1), 1).*x.^2,2);
    end
  end
  C(abs(C)<eps)=0;
  
end

function reccf = gpcf_linear_recappend(reccf, ri, gpcf)
%RECAPPEND Record append
%
%  Description
%    RECCF = GPCF_LINEAR_RECAPPEND(RECCF, RI, GPCF) takes a
%    covariance function record structure RECCF, record index RI
%    and covariance function structure GPCF with the current MCMC
%    samples of the parameters. Returns RECCF which contains all
%    the old samples and the current samples from GPCF. This 
%    subfunction is needed when using MCMC sampling (gp_mc).
%
%  See also
%    GP_MC and GP_MC -> RECAPPEND

  if nargin == 2
    % Initialize the record
    reccf.type = 'gpcf_linear';

    % Initialize parameters
    reccf.coeffSigma2= [];

    % Set the function handles
    reccf.fh.pak = @gpcf_linear_pak;
    reccf.fh.unpak = @gpcf_linear_unpak;
    reccf.fh.lp = @gpcf_linear_lp;
    reccf.fh.lpg = @gpcf_linear_lpg;
    reccf.fh.cfg = @gpcf_linear_cfg;
    reccf.fh.cov = @gpcf_linear_cov;
    reccf.fh.trcov  = @gpcf_linear_trcov;
    reccf.fh.trvar  = @gpcf_linear_trvar;
    reccf.fh.recappend = @gpcf_linear_recappend;
    reccf.p=[];
    reccf.p.coeffSigma2=[];
    if ~isempty(ri.p.coeffSigma2)
      reccf.p.coeffSigma2 = ri.p.coeffSigma2;
    end

  else
    % Append to the record
    gpp = gpcf.p;
    
    % record coeffSigma2
    reccf.coeffSigma2(ri,:)=gpcf.coeffSigma2;
    if isfield(gpp,'coeffSigma2') && ~isempty(gpp.coeffSigma2)
      reccf.p.coeffSigma2 = gpp.coeffSigma2.fh.recappend(reccf.p.coeffSigma2, ri, gpcf.p.coeffSigma2);
    end
  
    if isfield(gpcf, 'selectedVariables')
      reccf.selectedVariables = gpcf.selectedVariables;
    end
  end
end
