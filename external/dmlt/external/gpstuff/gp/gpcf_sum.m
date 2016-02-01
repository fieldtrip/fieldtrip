function gpcf = gpcf_sum(varargin)
%GPCF_SUM  Create a sum form covariance function
%
%  Description
%    GPCF = GPCF_SUM('cf', {GPCF_1, GPCF_2, ...}) 
%    creates a sum form covariance function
%          GPCF = GPCF_1 + GPCF_2 + ... + GPCF_N
%
%  See also
%    GP_SET, GPCF_*
  
% Copyright (c) 2009-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPCF_SUM';
  ip.addOptional('gpcf', [], @isstruct);
  ip.addParamValue('cf',[], @iscell);
  ip.parse(varargin{:});
  gpcf=ip.Results.gpcf;

  if isempty(gpcf)
    init=true;
    gpcf.type = 'gpcf_sum';
  else
    if ~isfield(gpcf,'type') && ~isequal(gpcf.type,'gpcf_sum')
      error('First argument does not seem to be a valid covariance function structure')
    end
    init=false;
  end
  
  if init || ~ismember('cf',ip.UsingDefaults)
    % Initialize parameters
    gpcf.cf = {};
    cfs=ip.Results.cf;
    if ~isempty(cfs)
      for i = 1:length(cfs)
        gpcf.cf{i} = cfs{i};
      end
    else
      error('At least one covariance function has to be given in cf');
    end
  end
  
  if init
    % Set the function handles to the subfunctions
    gpcf.fh.pak = @gpcf_sum_pak;
    gpcf.fh.unpak = @gpcf_sum_unpak;
    gpcf.fh.lp = @gpcf_sum_lp;
    gpcf.fh.lpg = @gpcf_sum_lpg;
    gpcf.fh.cfg = @gpcf_sum_cfg;
    gpcf.fh.ginput = @gpcf_sum_ginput;
    gpcf.fh.cov = @gpcf_sum_cov;
    gpcf.fh.trcov  = @gpcf_sum_trcov;
    gpcf.fh.trvar  = @gpcf_sum_trvar;
    gpcf.fh.recappend = @gpcf_sum_recappend;
  end

end

function [w, s] = gpcf_sum_pak(gpcf)
%GPCF_SUM_PAK  Combine GP covariance function parameters into one vector
%
%  Description
%    W = GPCF_SUM_PAK(GPCF, W) loops through all the covariance
%    functions and packs their parameters into one vector as
%    described in the respective functions. This is a mandatory 
%    subfunction used for example in energy and gradient computations.

%
%  See also
%    GPCF_SUM_UNPAK
  
  ncf = length(gpcf.cf);
  w = []; s = {};
  
  for i=1:ncf
    cf = gpcf.cf{i};
    [wi si] = cf.fh.pak(cf);
    w = [w wi];
    s = [s; si];
  end
end

function [gpcf, w] = gpcf_sum_unpak(gpcf, w)
%GPCF_SUM_UNPAK  Sets the covariance function parameters into
%                 the structures
%
%  Description
%    [GPCF, W] = GPCF_SUM_UNPAK(GPCF, W) loops through all the
%    covariance functions and unpacks their parameters from W to
%    each covariance function structure. This is a mandatory 
%    subfunction used for example in energy and gradient computations.

%
%  See also
%    GPCF_SUM_PAK
%
  ncf = length(gpcf.cf);
  
  for i=1:ncf
    cf = gpcf.cf{i};
    [cf, w] = cf.fh.unpak(cf, w);
    gpcf.cf{i} = cf;
  end

end

function lp = gpcf_sum_lp(gpcf)
%GPCF_SUM_LP  Evaluate the log prior of covariance function parameters
%
%  Description
%    LP = GPCF_SUM_LP(GPCF, X, T) takes a covariance function
%    structure GPCF and returns log(p(th)), where th collects the
%    parameters. This is a mandatory subfunction used for example 
%    in energy computations.
%
%  See also
%    GPCF_SUM_PAK, GPCF_SUM_UNPAK, GPCF_SUM_LPG, GP_E
  
  lp = 0;
  ncf = length(gpcf.cf);
  for i=1:ncf
    cf = gpcf.cf{i};
    lp = lp + cf.fh.lp(cf);
  end
  
end

function lpg = gpcf_sum_lpg(gpcf)
%GPCF_SUM_LPG  Evaluate gradient of the log prior with respect
%               to the parameters.
%
%  Description
%    LPG = GPCF_SUM_LPG(GPCF) takes a covariance function
%    structure GPCF and returns LPG = d log (p(th))/dth, where th
%    is the vector of parameters. This is a mandatory subfunction 
%    used for example in gradient computations.

%
%  See also
%    GPCF_SUM_PAK, GPCF_SUM_UNPAK, GPCF_SUM_LP, GP_G
  lpg = [];
  ncf = length(gpcf.cf);
  
  % Evaluate the gradients
  for i=1:ncf
    cf = gpcf.cf{i};
    lpg=[lpg cf.fh.lpg(cf)];
  end

end

function DKff = gpcf_sum_cfg(gpcf, x, x2, mask, i1)
%GPCF_SUM_CFG  Evaluate gradient of covariance function
%               with respect to the parameters.
%
%  Description
%    DKff = GPCF_SUM_CFG(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to th (cell array with matrix elements). This is a
%    mandatory subfunction used in gradient computations.
%
%    DKff = GPCF_SUM_CFG(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_SUM_CFG(GPCF, X, [], MASK) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the diagonal of gradients of covariance matrix
%    Kff = k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_SUM_CFG(GPCF, X, X2, [], i) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2), or k(X,X) if X2 is empty, with respect to ith 
%    hyperparameter. This subfunction is needed when using 
%    memory save option in gp_set.
%
%  See also
%    GPCF_SUM_PAK, GPCF_SUM_UNPAK, GPCF_SUM_LP, GP_G

  [n, m] =size(x);
  ncf = length(gpcf.cf);

  DKff = {};
 
  if nargin==5
    % Use memory save option
    savememory=1;
    i3=0;
    for k=1:ncf
      % Number of hyperparameters for each covariance funtion
      cf=gpcf.cf{k};
      i3(k)=cf.fh.cfg(cf,[],[],[],0);
    end
    if i1==0
      % Return number of hyperparameters
      DKff=sum(i3);
      return
    end
    % Help indices
    i3=cumsum(i3);
    ind=find(cumsum(i3 >= i1)==1);
    if ind>1
      i1=[ind i1-i3(ind-1)];
    else
      i1=[ind i1];
    end
    i2=i1(1);
  else
    savememory=0;
    i2=1:ncf;
  end

  % Evaluate: DKff{1} = d Kff / d magnSigma2
  %           DKff{2} = d Kff / d lengthScale
  % NOTE! Here we have already taken into account that the parameters are transformed
  % through log() and thus dK/dlog(p) = p * dK/dp

  % evaluate the gradient for training covariance
  if nargin == 2 || (isempty(x2) && isempty(mask))
        
%     % evaluate the individual covariance functions
%     for i=1:ncf
%       cf = gpcf.cf{i};
%       C{i} = cf.fh.trcov(cf, x);
%     end
    
    % Evaluate the gradients
    ind = 1:ncf;
    DKff = {};
    for i=i2
      cf = gpcf.cf{i};
      if ~savememory
        DK = cf.fh.cfg(cf, x);
      else
        DK{1} = cf.fh.cfg(cf,x,[],[],i1(2));
      end
      DKff = [DKff DK];

    end
    
    % Evaluate the gradient of non-symmetric covariance (e.g. K_fu)
  elseif nargin == 3 || isempty(mask)
    if size(x,2) ~= size(x2,2)
      error('gpcf_sum -> _ghyper: The number of columns in x and x2 has to be the same. ')
    end
        
    % Evaluate the gradients
    ind = 1:ncf;
    DKff = {};
    for i=i2
      cf = gpcf.cf{i};
      if ~savememory
        DK = cf.fh.cfg(cf, x, x2);
      else
        DK{1} = cf.fh.cfg(cf,x,x2,[],i1(2));
      end
      DKff = [DKff DK];

    end

    
    
    % Evaluate: DKff{1}    = d mask(Kff,I) / d magnSigma2
    %           DKff{2...} = d mask(Kff,I) / d lengthScale
  elseif nargin == 4 || nargin == 5
    
    % Evaluate the gradients
    ind = 1:ncf;
    DKff = {};
    for i=i2
      cf = gpcf.cf{i};
      if savememory
        DK = cf.fh.cfg(cf, x, [], 1, i1(2));
      else
        DK = cf.fh.cfg(cf, x, [], 1);
      end
      DKff = [DKff DK];

    end
  end
  if savememory
    DKff=DKff{1};
  end
end


function DKff = gpcf_sum_ginput(gpcf, x, x2, i1)
%GPCF_SUM_GINPUT  Evaluate gradient of covariance function with 
%                  respect to x
%
%  Description
%    DKff = GPCF_SUM_GINPUT(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to X (cell array with matrix elements). This subfunction 
%    is needed when computing gradients with respect to inducing 
%    inputs in sparse approximations.
%
%    DKff = GPCF_SUM_GINPUT(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_SUM_GINPUT(GPCF, X, X2, i) takes a covariance
%    function structure GPCF, a matrix X of input vectors
%    and returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2), or k(X,X) if X2 is empty, with respect to ith 
%    covariate in X. This subfunction is needed when using 
%    memory save option in gp_set.
%
%  See also
%    GPCF_SUM_PAK, GPCF_SUM_UNPAK, GPCF_SUM_LP, GP_G
  
  [n, m] =size(x);
  
  ncf = length(gpcf.cf);
  
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
  if ~savememory
    DKff=cellfun(@(a) zeros(n,n), cell(1,numel(x)), 'UniformOutput', 0);
  else
    DKff=cellfun(@(a) zeros(n,n), cell(1,n), 'UniformOutput', 0);
  end
  if nargin == 2 || isempty(x2)
    
    % Evaluate the gradients
    ind = 1:ncf;
    for i=1:ncf
      cf = gpcf.cf{i};
      if ~savememory
        DK = cf.fh.ginput(cf, x);
      else
        DK = cf.fh.ginput(cf,x,[],i1);
      end
      for j=1:length(DK)
        DKff{j} = DKff{j} + DK{j};
      end

    end

    % Evaluate the gradient of non-symmetric covariance (e.g. K_fu)
  elseif nargin == 3 || nargin == 4
    if size(x,2) ~= size(x2,2)
      error('gpcf_sum -> _ghyper: The number of columns in x and x2 has to be the same. ')
    end  
    
    % Evaluate the gradients
    for i=1:ncf
      cf = gpcf.cf{i};
      if ~savememory
        DK = cf.fh.ginput(cf, x, x2);
      else
        DK = cf.fh.ginput(cf,x,x2,i1);
      end
      for j=1:length(DK)
        DKff{j} = DKff{j} + DK{j};
      end
      
    end
  end
  
end


function C = gpcf_sum_cov(gpcf, x1, x2)
%GP_SUM_COV  Evaluate covariance matrix between two input vectors
%
%  Description         
%    C = GP_SUM_COV(GP, TX, X) takes in covariance function of a
%    Gaussian process GP and two matrixes TX and X that contain
%    input vectors to GP. Returns covariance matrix C. Every
%    element ij of C contains covariance between inputs i in TX
%    and j in X. This is a mandatory subfunction used for example in
%    prediction and energy computations.
%
%
%  See also
%    GPCF_SUM_TRCOV, GPCF_SUM_TRVAR, GP_COV, GP_TRCOV
  
  if isempty(x2)
    x2=x1;
  end
  [n1,m1]=size(x1);
  [n2,m2]=size(x2);

  if m1~=m2
    error('the number of columns of X1 and X2 has to be same')
  end

  ncf = length(gpcf.cf);
  
  % evaluate the individual covariance functions
  C = 0;
  for i=1:ncf
    cf = gpcf.cf{i};
    C = C + cf.fh.cov(cf, x1, x2);
  end        
end

function C = gpcf_sum_trcov(gpcf, x)
%GP_SUM_TRCOV     Evaluate training covariance matrix of inputs
%
%  Description
%    C = GP_SUM_TRCOV(GP, TX) takes in covariance function of a
%    Gaussian process GP and matrix TX that contains training
%    input vectors. Returns covariance matrix C. Every element ij
%    of C contains covariance between inputs i and j in TX. This 
%    is a mandatory subfunction used for example in prediction and 
%    energy computations.
%
%  See also
%    GPCF_SUM_COV, GPCF_SUM_TRVAR, GP_COV, GP_TRCOV
  ncf = length(gpcf.cf);
  
  % evaluate the individual covariance functions
  C = 0;
  for i=1:ncf
    cf = gpcf.cf{i};
    C = C + cf.fh.trcov(cf, x);
  end
end

function C = gpcf_sum_trvar(gpcf, x)
% GP_SUM_TRVAR     Evaluate training variance vector
%
%  Description
%    C = GP_SUM_TRVAR(GPCF, TX) takes in covariance function of
%    a Gaussian process GPCF and matrix TX that contains training
%    inputs. Returns variance vector C. Every element i of C
%    contains variance of input i in TX. This is a mandatory 
%    subfunction used for example in prediction and energy computations.
%
%  See also
%    GPCF_SUM_COV, GP_COV, GP_TRCOV


  ncf = length(gpcf.cf);
  
  % evaluate the individual covariance functions
  C = 0;
  for i=1:ncf
    cf = gpcf.cf{i};
    C = C + cf.fh.trvar(cf, x);
  end
end

function reccf = gpcf_sum_recappend(reccf, ri, gpcf)
%RECAPPEND  Record append
%
%  Description
%    RECCF = GPCF_SUM_RECAPPEND(RECCF, RI, GPCF) takes a
%    covariance function record structure RECCF, record index RI
%    and covariance function structure GPCF with the current MCMC
%    samples of the parameters. Returns RECCF which contains all
%    the old samples and the current samples from GPCF. This 
%    subfunction is needed when using MCMC sampling (gp_mc).
%
%  See also
%    GP_MC, GP_MC->RECAPPEND
  
  if nargin == 2
    % Initialize the record
    reccf.type = 'gpcf_sum';

    % Initialize parameters
    ncf = length(ri.cf);
    for i=1:ncf
      cf = ri.cf{i};
      reccf.cf{i} = cf.fh.recappend([], ri.cf{i});
    end
    
    % Set the function handles
    reccf.fh.pak = @gpcf_sum_pak;
    reccf.fh.unpak = @gpcf_sum_unpak;
    reccf.fh.e = @gpcf_sum_lp;
    reccf.fh.lpg = @gpcf_sum_lpg;
    reccf.fh.cfg = @gpcf_sum_cfg;
    reccf.fh.cov = @gpcf_sum_cov;
    reccf.fh.trcov  = @gpcf_sum_trcov;
    reccf.fh.trvar  = @gpcf_sum_trvar;
    reccf.fh.recappend = @gpcf_sum_recappend;
  else
    % Append to the record
    
    % Loop over all of the covariance functions
    ncf = length(gpcf.cf);
    for i=1:ncf
      cf = gpcf.cf{i};
      reccf.cf{i} = cf.fh.recappend(reccf.cf{i}, ri, cf);
    end
  end
end

