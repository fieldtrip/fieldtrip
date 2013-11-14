function gpcf = gpcf_neuralnetwork(varargin)
%GPCF_NEURALNETWORK  Create a neural network covariance function
%
%  Description
%    GPCF = GPCF_NEURALNETWORK('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates neural network covariance function structure in which
%    the named parameters have the specified values. Any
%    unspecified parameters are set to default values.
%
%    GPCF = GPCF_NEURALNETWORK(GPCF,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a covariance function structure with the named
%    parameters altered with the specified values.
%  
%    Parameters for neural network covariance function [default]
%      biasSigma2         - prior variance for bias in neural network [0.1]
%      weightSigma2       - prior variance for weights in neural network [10]
%                           This can be either scalar corresponding
%                           to a common prior variance or vector
%                           defining own prior variance for each
%                           input.
%      biasSigma2_prior   - prior structure for magnSigma2 [prior_logunif]
%      weightSigma2_prior - prior structure for weightSigma2 [prior_logunif]
%      selectedVariables  - vector defining which inputs are used
%
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc.
%
%  See also
%    GP_SET, GPCF_*, PRIOR_*
  
% Copyright (c) 2007-2009 Jarno Vanhatalo
% Copyright (c) 2009 Jaakko Riihimaki
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPCF_NEURALNETWORK';
  ip.addOptional('gpcf', [], @isstruct);
  ip.addParamValue('biasSigma2',0.1, @(x) isscalar(x) && x>0);
  ip.addParamValue('weightSigma2',10, @(x) isvector(x) && all(x>0));
  ip.addParamValue('biasSigma2_prior',prior_logunif, @(x) isstruct(x) || isempty(x));
  ip.addParamValue('weightSigma2_prior',prior_logunif, @(x) isstruct(x) || isempty(x));
  ip.addParamValue('selectedVariables',[], @(x) isvector(x) && all(x>0));
  ip.parse(varargin{:});
  gpcf=ip.Results.gpcf;
  
  if isempty(gpcf)
    init=true;
    gpcf.type = 'gpcf_neuralnetwork';
  else
    if ~isfield(gpcf,'type') && ~isequal(gpcf.type,'gpcf_neuralnetwork')
      error('First argument does not seem to be a valid covariance function structure')
    end
    init=false;
  end
  
  % Initialize parameters
  if init || ~ismember('biasSigma2',ip.UsingDefaults)
    gpcf.biasSigma2=ip.Results.biasSigma2;
  end
  if init || ~ismember('weightSigma2',ip.UsingDefaults)
    gpcf.weightSigma2=ip.Results.weightSigma2;
  end

  % Initialize prior structure
  if init
    gpcf.p=[];
  end
  if init || ~ismember('biasSigma2_prior',ip.UsingDefaults)
    gpcf.p.biasSigma2=ip.Results.biasSigma2_prior;
  end
  if init || ~ismember('weightSigma2_prior',ip.UsingDefaults)
    gpcf.p.weightSigma2=ip.Results.weightSigma2_prior;
  end
  if ~ismember('selectedVariables',ip.UsingDefaults)
    selectedVariables=ip.Results.selectedVariables;
    if ~isempty(selectedVariables)
      gpcf.selectedVariables = selectedVariables;
    end
  end
  
  if init
    % Set the function handles to the subfunctions
    gpcf.fh.pak = @gpcf_neuralnetwork_pak;
    gpcf.fh.unpak = @gpcf_neuralnetwork_unpak;
    gpcf.fh.lp = @gpcf_neuralnetwork_lp;
    gpcf.fh.lpg = @gpcf_neuralnetwork_lpg;
    gpcf.fh.cfg = @gpcf_neuralnetwork_cfg;
    gpcf.fh.ginput = @gpcf_neuralnetwork_ginput;
    gpcf.fh.cov = @gpcf_neuralnetwork_cov;
    gpcf.fh.trcov  = @gpcf_neuralnetwork_trcov;
    gpcf.fh.trvar  = @gpcf_neuralnetwork_trvar;
    gpcf.fh.recappend = @gpcf_neuralnetwork_recappend;
  end

end

function [w, s] = gpcf_neuralnetwork_pak(gpcf, w)
%GPCF_NEURALNETWORK_PAK  Combine GP covariance function parameters
%                        into one vector
%
%  Description
%   W = GPCF_NEURALNETWORK_PAK(GPCF) takes a covariance function
%   structure GPCF and combines the covariance function parameters
%   and their hyperparameters into a single row vector W. This is a 
%   mandatory subfunction used for example in energy and gradient 
%   computations.
%
%       w = [ log(gpcf.biasSigma2)
%             (hyperparameters of gpcf.biasSigma2) 
%             log(gpcf.weightSigma2(:))
%             (hyperparameters of gpcf.weightSigma2)]'
%     
%
%  See also
%   GPCF_NEURALNETWORK_UNPAK

  i1=0;i2=1;
  w = []; s = {};
  
  if ~isempty(gpcf.p.biasSigma2)
    w = [w log(gpcf.biasSigma2)];
    s = [s; 'log(neuralnetwork.biasSigma2)'];
    
    % Hyperparameters of magnSigma2
    [wh sh] = gpcf.p.biasSigma2.fh.pak(gpcf.p.biasSigma2);
    w = [w wh];
    s = [s; sh];
  end        
  
  if ~isempty(gpcf.p.weightSigma2)
    w = [w log(gpcf.weightSigma2)];
    if numel(gpcf.weightSigma2)>1
      s = [s; sprintf('log(neuralnetwork.weightSigma2 x %d)',numel(gpcf.weightSigma2))];
    else
      s = [s; 'log(neuralnetwork.weightSigma2)'];
    end
    
    % Hyperparameters of lengthScale
    [wh sh] = gpcf.p.weightSigma2.fh.pak(gpcf.p.weightSigma2);
    w = [w wh];
    s = [s; sh];
  end

end


function [gpcf, w] = gpcf_neuralnetwork_unpak(gpcf, w)
%GPCF_NEURALNETWORK_UNPAK  Sets the covariance function parameters 
%                          into the structure
%
%  Description
%    [GPCF, W] = GPCF_NEURALNETWORK_UNPAK(GPCF, W) takes a
%    covariance function structure GPCF and a hyper-parameter
%    vector W, and returns a covariance function structure
%    identical to the input, except that the covariance
%    hyper-parameters have been set to the values in W. Deletes
%    the values set to GPCF from W and returns the modified W.
%    This is a mandatory subfunction used for example in energy 
%    and gradient computations.
%
%    Assignment is inverse of  
%       w = [ log(gpcf.coeffSigma2)
%             (hyperparameters of gpcf.coeffSigma2)]'
%
%  See also
%   GPCF_NEURALNETWORK_PAK
  
  gpp=gpcf.p;
  if ~isempty(gpp.biasSigma2)
    i1=1;
    gpcf.biasSigma2 = exp(w(i1));
    w = w(i1+1:end);
  end

  if ~isempty(gpp.weightSigma2)
    i2=length(gpcf.weightSigma2);
    i1=1;
    gpcf.weightSigma2 = exp(w(i1:i2));
    w = w(i2+1:end);
    
    % Hyperparameters of lengthScale
    [p, w] = gpcf.p.weightSigma2.fh.unpak(gpcf.p.weightSigma2, w);
    gpcf.p.weightSigma2 = p;
  end
  
  if ~isempty(gpp.biasSigma2)
    % Hyperparameters of magnSigma2
    [p, w] = gpcf.p.biasSigma2.fh.unpak(gpcf.p.biasSigma2, w);
    gpcf.p.biasSigma2 = p;
  end
end

function lp = gpcf_neuralnetwork_lp(gpcf)
%GPCF_NEURALNETWORK_LP  Evaluate the log prior of covariance
%                       function parameters
%
%  Description
%    LP = GPCF_NEURALNETWORK_LP(GPCF) takes a covariance function
%    structure GPCF and returns log(p(th)), where th collects the
%    parameters. This is a mandatory subfunction used for example 
%    in energy computations.
%
%  See also
%   GPCF_NEURALNETWORK_PAK, GPCF_NEURALNETWORK_UNPAK,
%   GPCF_NEURALNETWORK_LPG, GP_E

  lp = 0;
  gpp=gpcf.p;

  if ~isempty(gpp.biasSigma2)
    lp = gpp.biasSigma2.fh.lp(gpcf.biasSigma2, gpp.biasSigma2) +log(gpcf.biasSigma2);
  end
  if ~isempty(gpp.weightSigma2)
    lp = lp +gpp.weightSigma2.fh.lp(gpcf.weightSigma2, gpp.weightSigma2) +sum(log(gpcf.weightSigma2));
  end

end

function lpg = gpcf_neuralnetwork_lpg(gpcf)
%GPCF_NEURALNETWORK_LPG  Evaluate gradient of the log prior with respect
%                 to the parameters.
%
%  Description
%    LPG = GPCF_NEURALNETWORK_LPG(GPCF) takes a covariance function
%    structure GPCF and returns LPG = d log (p(th))/dth, where th
%    is the vector of parameters. This is a mandatory subfunction 
%    used for example in gradient computations.
%
%  See also
%    GPCF_NEURALNETWORK_PAK, GPCF_NEURALNETWORK_UNPAK,
%    GPCF_NEURALNETWORK_LP, GP_G

  lpg = [];
  gpp=gpcf.p;
  
  if ~isempty(gpcf.p.biasSigma2)            
    lpgs = gpp.biasSigma2.fh.lpg(gpcf.biasSigma2, gpp.biasSigma2);
    lpg = [lpg lpgs(1).*gpcf.biasSigma2+1 lpgs(2:end)];
  end
  if ~isempty(gpcf.p.weightSigma2)
    lll = length(gpcf.weightSigma2);
    lpgs = gpp.weightSigma2.fh.lpg(gpcf.weightSigma2, gpp.weightSigma2);
    lpg = [lpg lpgs(1:lll).*gpcf.weightSigma2+1 lpgs(lll+1:end)];
  end
end

function DKff = gpcf_neuralnetwork_cfg(gpcf, x, x2, mask, i1)
%GPCF_NEURALNETWORK_CFG  Evaluate gradient of covariance function
%                        with respect to the parameters
%
%  Description
%    DKff = GPCF_NEURALNETWORK_CFG(GPCF, X) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X) with respect to th (cell array with matrix elements).
%    This is a mandatory subfunction used in gradient computations.
%
%    DKff = GPCF_NEURALNETWORK_CFG(GPCF, X, X2) takes a
%    covariance function structure GPCF, a matrix X of input
%    vectors and returns DKff, the gradients of covariance matrix
%    Kff = k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_NEURALNETWORK_CFG(GPCF, X, [], MASK) takes a
%    covariance function structure GPCF, a matrix X of input
%    vectors and returns DKff, the diagonal of gradients of
%    covariance matrix Kff = k(X,X2) with respect to th (cell
%    array with matrix elements). This subfunction is needed 
%    when using sparse approximations (e.g. FIC).
%
%    DKff = GPCF_NEURALNETWORK_CFG(GPCF,X,X2,[],i) takes a
%    covariance function structure GPCF, a matrix X of input
%    vectors and returns DKff, the gradients of covariance matrix
%    Kff = k(X,X2) with respect to ith hyperparameter(matrix).
%    5th input parameter can also be used without X2. If i = 0, 
%    number of hyperparameters is returned. This subfunction is
%    needed when using memory save option in gp_set.
%  
%
%  See also
%    GPCF_NEURALNETWORK_PAK, GPCF_NEURALNETWORK_UNPAK,
%    GPCF_NEURALNETWORK_LP, GP_G
  
  gpp=gpcf.p;
  
  if isfield(gpcf, 'selectedVariables') && ~isempty(x)
    x=x(:,gpcf.selectedVariables);
    if nargin == 3
      x2=x2(:,gpcf.selectedVariables);
    end
  end
  [n, m] =size(x);
    
  if nargin==5
    % Use memory save option
    if i1==0
      % Return number of hyperparameters
      if ~isempty(gpcf.p.biasSigma2)
        i=1;
      end
      if ~isempty(gpcf.p.weightSigma2)
        i=i+length(gpcf.weightSigma2);
      end
      DKff=i;
      return
    end
    savememory=1;
  else
    savememory=0;
    i1=1:m;
  end  
  
  DKff = {};
  gprior = [];
  
  % Evaluate: DKff{1} = d Kff / d biasSigma2
  %           DKff{2} = d Kff / d weightSigma2
  % NOTE! Here we have already taken into account that the parameters
  % are transformed through log() and thus dK/dlog(p) = p * dK/dp
  
  % evaluate the gradient for training covariance
  if nargin == 2 || (isempty(x2) && isempty(mask))
    
    x_aug=[ones(size(x,1),1) x];
    
    if length(gpcf.weightSigma2) == 1
      % In the case of an isotropic NEURALNETWORK
      s = gpcf.weightSigma2*ones(1,m);
    else
      s = gpcf.weightSigma2;
    end
    
    S_nom=2*x_aug*diag([gpcf.biasSigma2 s])*x_aug';
    
    S_den_tmp=(2*sum(repmat([gpcf.biasSigma2 s], n, 1).*x_aug.^2,2)+1);
    S_den2=S_den_tmp*S_den_tmp';
    S_den=sqrt(S_den2);
    
    C_tmp=2/pi./sqrt(1-(S_nom./S_den).^2);
    % C(abs(C)<=eps) = 0;
    C_tmp = (C_tmp+C_tmp')./2;
    ii1 = 0;    
    
    if ~savememory || i1==1
      bnom_g=2*ones(n);
      bden_g=(0.5./S_den).*(bnom_g.*repmat(S_den_tmp',n,1)+repmat(S_den_tmp,1,n).*bnom_g);
      bg=gpcf.biasSigma2*C_tmp.*(bnom_g.*S_den-bden_g.*S_nom)./S_den2;
      
      if ~isempty(gpcf.p.biasSigma2)
        ii1 = ii1+1;
        DKff{ii1}=(bg+bg')/2;
      end
      if savememory
        DKff=DKff{ii1};
        return
      end
    elseif savememory
      i1=i1-1;
    end
    
    if ~isempty(gpcf.p.weightSigma2)
      if length(gpcf.weightSigma2) == 1
        wnom_g=2*x*x';
        tmp_g=sum(2*x.^2,2);
        wden_g=0.5./S_den.*(tmp_g*S_den_tmp'+S_den_tmp*tmp_g');
        wg=s(1)*C_tmp.*(wnom_g.*S_den-wden_g.*S_nom)./S_den2;
        
        ii1 = ii1+1;
        DKff{ii1}=(wg+wg')/2;
      else
        for d1=i1
          wnom_g=2*x(:,d1)*x(:,d1)';
          tmp_g=2*x(:,d1).^2;
          tmp=tmp_g*S_den_tmp';
          wden_g=0.5./S_den.*(tmp+tmp');
          wg=s(d1)*C_tmp.*(wnom_g.*S_den-wden_g.*S_nom)./S_den2;
          
          ii1 = ii1+1;
          DKff{ii1}=(wg+wg')/2;
        end
      end
    end
    
    % Evaluate the gradient of non-symmetric covariance (e.g. K_fu)
  elseif nargin == 3 || isempty(mask)
    
    if size(x,2) ~= size(x2,2)
      error('gpcf_neuralnetwork -> _ghyper: The number of columns in x and x2 has to be the same. ')
    end
    
    n2 =size(x2,1);
    
    x_aug=[ones(size(x,1),1) x];
    x_aug2=[ones(size(x2,1),1) x2];
    
    if length(gpcf.weightSigma2) == 1
      % In the case of an isotropic NEURALNETWORK
      s = gpcf.weightSigma2*ones(1,m);
    else
      s = gpcf.weightSigma2;
    end
    
    S_nom=2*x_aug*diag([gpcf.biasSigma2 s])*x_aug2';
    
    S_den_tmp1=(2*sum(repmat([gpcf.biasSigma2 s], n, 1).*x_aug.^2,2)+1);
    S_den_tmp2=(2*sum(repmat([gpcf.biasSigma2 s], n2, 1).*x_aug2.^2,2)+1);
    
    S_den2=S_den_tmp1*S_den_tmp2';
    S_den=sqrt(S_den2);
    
    C_tmp=2/pi./sqrt(1-(S_nom./S_den).^2);
    %C(abs(C)<=eps) = 0;
    ii1 = 0;

    if ~savememory || i1==1
      bnom_g=2*ones(n, n2);
      bden_g=(0.5./S_den).*(bnom_g.*repmat(S_den_tmp2',n,1)+repmat(S_den_tmp1,1,n2).*bnom_g);
      
      if ~isempty(gpcf.p.biasSigma2)
        ii1 = ii1 + 1;
        DKff{ii1}=gpcf.biasSigma2*C_tmp.*(bnom_g.*S_den-bden_g.*S_nom)./S_den2;
      end
      if savememory
        DKff=DKff{ii1};
        return
      end
    elseif savememory
      i1=i1-1;
    end
    
    if ~isempty(gpcf.p.weightSigma2)
      if length(gpcf.weightSigma2) == 1
        wnom_g=2*x*x2';
        tmp_g1=sum(2*x.^2,2);
        tmp_g2=sum(2*x2.^2,2);
        wden_g=0.5./S_den.*(tmp_g1*S_den_tmp2'+S_den_tmp1*tmp_g2');
        
        ii1 = ii1 + 1;
        DKff{ii1}=s(1)*C_tmp.*(wnom_g.*S_den-wden_g.*S_nom)./S_den2;
      else
        for d1=i1
          wnom_g=2*x(:,d1)*x2(:,d1)';
          tmp_g1=2*x(:,d1).^2;
          tmp_g2=2*x2(:,d1).^2;
          wden_g=0.5./S_den.*(tmp_g1*S_den_tmp2'+S_den_tmp1*tmp_g2');
          
          ii1 = ii1 + 1;
          DKff{ii1}=s(d1)*C_tmp.*(wnom_g.*S_den-wden_g.*S_nom)./S_den2;
        end
      end
    end

    % Evaluate: DKff{1}    = d mask(Kff,I) / d biasSigma2
    %           DKff{2...} = d mask(Kff,I) / d weightSigma2
  elseif nargin == 4 || nargin == 5
    
    x_aug=[ones(size(x,1),1) x];
    
    if length(gpcf.weightSigma2) == 1
      % In the case of an isotropic NEURALNETWORK
      s = gpcf.weightSigma2*ones(1,m);
    else
      s = gpcf.weightSigma2;
    end
    
    S_nom=2*sum(repmat([gpcf.biasSigma2 s],n,1).*x_aug.^2,2);
    
    S_den=(S_nom+1);
    S_den2=S_den.^2;
    
    C_tmp=2/pi./sqrt(1-(S_nom./S_den).^2);
    %C(abs(C)<=eps) = 0;
    
    bnom_g=2*ones(n,1);
    bden_g=(0.5./S_den).*(2*bnom_g.*S_den);
    
    ii1 = 0;
    if ~isempty(gpcf.p.biasSigma2) && (~savememory || all(i1==1))
      ii1 = ii1 + 1;
      DKff{ii1}=gpcf.biasSigma2*C_tmp.*(bnom_g.*S_den-bden_g.*S_nom)./S_den2;
    end
    if savememory
      if i1==1
        DKff=DKff{1};
        return
      end
      i1=i1-1;
    end
    
    if ~isempty(gpcf.p.weightSigma2)
      if length(gpcf.weightSigma2) == 1
        wnom_g=sum(2*x.^2,2);
        wden_g=0.5./S_den.*(2*wnom_g.*S_den);
        
        ii1 = ii1+1;
        DKff{ii1}=s(1)*C_tmp.*(wnom_g.*S_den-wden_g.*S_nom)./S_den2;
      else
        for d1=i1
          wnom_g=2*x(:,d1).^2;
          wden_g=0.5./S_den.*(2*wnom_g.*S_den);
          
          ii1 = ii1+1;                        
          DKff{ii1}=s(d1)*C_tmp.*(wnom_g.*S_den-wden_g.*S_nom)./S_den2;
        end
      end
    end
  end
  if savememory
    DKff=DKff{1};
  end
  
end


function DKff = gpcf_neuralnetwork_ginput(gpcf, x, x2, i1)
%GPCF_NEURALNETWORK_GINPUT  Evaluate gradient of covariance function with 
%                           respect to x.
%
%  Description
%    DKff = GPCF_NEURALNETWORK_GINPUT(GPCF, X) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_NEURALNETWORK_GINPUT(GPCF, X, X2) takes a
%    covariance function structure GPCF, a matrix X of input
%    vectors and returns DKff, the gradients of covariance matrix
%    Kff = k(X,X2) with respect to X (cell array with matrix
%    elements). This subfunction is needed when computing gradients 
%    with respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_NEURALNETWORK_GINPUT(GPCF, X, X2, i) takes a
%    covariance function structure GPCF, a matrix X of input
%    vectors and returns DKff, the gradients of covariance matrix
%    Kff = k(X,X2) with respect to ith covariate in  X (matrix).
%    This subfunction is needed when using memory save option in
%    gp_set.
%
%  See also
%    GPCF_NEURALNETWORK_PAK, GPCF_NEURALNETWORK_UNPAK,
%    GPCF_NEURALNETWORK_LP, GP_G
  
  if isfield(gpcf, 'selectedVariables')
    x=x(:,gpcf.selectedVariables); 
    if nargin == 3
      x2=x2(:,gpcf.selectedVariables); 
    end
  end
  
  [n, m] =size(x);
  
  if nargin==4
    % Use memory save option
    if i1==0
      % Return number of covariates
      DKff=m;
      return
    end
  else
    i1=1:m;
  end
  
  if nargin == 2 || isempty(x2) 
    
    if length(gpcf.weightSigma2) == 1
      % In the case of an isotropic NEURALNETWORK
      s = gpcf.weightSigma2*ones(1,m);
    else
      s = gpcf.weightSigma2;
    end
    
    x_aug=[ones(size(x,1),1) x];
    
    S_nom=2*x_aug*diag([gpcf.biasSigma2 s])*x_aug';
    S_den_tmp=(2*sum(repmat([gpcf.biasSigma2 s], n, 1).*x_aug.^2,2)+1);
    
    S_den2=S_den_tmp*S_den_tmp';
    S_den=sqrt(S_den2);
    
    C_tmp=2/pi./sqrt(1-(S_nom./S_den).^2);
    %C(abs(C)<=eps) = 0;
    C_tmp = (C_tmp+C_tmp')./2;
    
    ii1=0;
    for d1=i1
      for j=1:n
        
        DK = zeros(n);
        DK(j,:)=s(d1)*x(:,d1)';
        DK = DK + DK';
        inom_g=2*DK;
        
        tmp_g=zeros(n);
        tmp_g(j,:)=2*s(d1)*2*x(j,d1)*S_den_tmp';
        tmp_g=tmp_g+tmp_g';
        
        iden_g=0.5./S_den.*(tmp_g);
        
        ii1=ii1+1;
        DKff{ii1}=C_tmp.*(inom_g.*S_den-iden_g.*S_nom)./S_den2;
      end
    end
    
  elseif nargin == 3 || nargin == 4
    
    if length(gpcf.weightSigma2) == 1
      % In the case of an isotropic NEURALNETWORK
      s = gpcf.weightSigma2*ones(1,m);
    else
      s = gpcf.weightSigma2;
    end
    
    n2 =size(x2,1);
    
    x_aug=[ones(size(x,1),1) x];
    x_aug2=[ones(size(x2,1),1) x2];
    
    S_nom=2*x_aug*diag([gpcf.biasSigma2 s])*x_aug2';
    
    S_den_tmp1=(2*sum(repmat([gpcf.biasSigma2 s], n, 1).*x_aug.^2,2)+1);
    S_den_tmp2=(2*sum(repmat([gpcf.biasSigma2 s], n2, 1).*x_aug2.^2,2)+1);
    
    S_den2=S_den_tmp1*S_den_tmp2';
    S_den=sqrt(S_den2);
    
    C_tmp=2/pi./sqrt(1-(S_nom./S_den).^2);
    % C(abs(C)<=eps) = 0;
    
    ii1 = 0;
    for d1=i1
      for j = 1:n
        
        DK = zeros(n, n2);
        DK(j,:)=s(d1)*x2(:,d1)';
        inom_g=2*DK;
        
        tmp_g=zeros(n, n2);
        tmp_g(j,:)=2*s(d1)*2*x(j,d1)*S_den_tmp2';
        
        iden_g=0.5./S_den.*(tmp_g);
        
        ii1=ii1+1;
        DKff{ii1}=C_tmp.*(inom_g.*S_den-iden_g.*S_nom)./S_den2;
      end
    end
  end
end


function C = gpcf_neuralnetwork_cov(gpcf, x1, x2, varargin)
%GP_NEURALNETWORK_COV  Evaluate covariance matrix between two input vectors
%
%  Description         
%    C = GP_NEURALNETWORK_COV(GP, TX, X) takes in covariance
%    function of a Gaussian process GP and two matrixes TX and X
%    that contain input vectors to GP. Returns covariance matrix
%    C. Every element ij of C contains covariance between inputs
%    i in TX and j in X. This is a mandatory subfunction used for 
%    example in prediction and energy computations.
%
%
%  See also
%    GPCF_NEURALNETWORK_TRCOV, GPCF_NEURALNETWORK_TRVAR, GP_COV,
%    GP_TRCOV
  
  if isfield(gpcf, 'selectedVariables')
    x1=x1(:,gpcf.selectedVariables); 
    if nargin == 3
      x2=x2(:,gpcf.selectedVariables); 
    end
  end
  
  
  if isempty(x2)
    x2=x1;
  end
  
  [n1,m1]=size(x1);
  [n2,m2]=size(x2);
  
  if m1~=m2
    error('the number of columns of X1 and X2 has to be same')
  end
  
  x_aug1=[ones(n1,1) x1];
  x_aug2=[ones(n2,1) x2];
  
  if length(gpcf.weightSigma2) == 1
    % In the case of an isotropic NEURALNETWORK
    s = gpcf.weightSigma2*ones(1,m1);
  else
    s = gpcf.weightSigma2;
  end
  
  S_nom=2*x_aug1*diag([gpcf.biasSigma2 s])*x_aug2';
  
  S_den_tmp1=(2*sum(repmat([gpcf.biasSigma2 s], n1, 1).*x_aug1.^2,2)+1);
  S_den_tmp2=(2*sum(repmat([gpcf.biasSigma2 s], n2, 1).*x_aug2.^2,2)+1);
  S_den2=S_den_tmp1*S_den_tmp2';
  
  C=2/pi*asin(S_nom./sqrt(S_den2));
  
  C(abs(C)<=eps) = 0;
end


function C = gpcf_neuralnetwork_trcov(gpcf, x)
%GP_NEURALNETWORK_TRCOV  Evaluate training covariance matrix of inputs
%
%  Description
%    C = GP_NEURALNETWORK_TRCOV(GP, TX) takes in covariance
%    function of a Gaussian process GP and matrix TX that
%    contains training input vectors. Returns covariance matrix
%    C. Every element ij of C contains covariance between inputs
%    i and j in TX. This is a mandatory subfunction used for 
%    example in prediction and energy computations.
%
%  See also
%    GPCF_NEURALNETWORK_COV, GPCF_NEURALNETWORK_TRVAR, GP_COV,
%    GP_TRCOV
  
  if isfield(gpcf, 'selectedVariables')
    x=x(:,gpcf.selectedVariables); 
  end
  
  [n,m]=size(x);
  x_aug=[ones(n,1) x];
  
  if length(gpcf.weightSigma2) == 1
    % In the case of an isotropic NEURALNETWORK
    s = gpcf.weightSigma2*ones(1,m);
  else
    s = gpcf.weightSigma2;
  end
  
  S_nom=2*x_aug*diag([gpcf.biasSigma2 s])*x_aug';
  
  S_den_tmp=(2*sum(repmat([gpcf.biasSigma2 s], n, 1).*x_aug.^2,2)+1);        
  S_den2=S_den_tmp*S_den_tmp';
  
  C=2/pi*asin(S_nom./sqrt(S_den2));
  
  C(abs(C)<=eps) = 0;
  C = (C+C')./2;
  
end

function C = gpcf_neuralnetwork_trvar(gpcf, x)
%GP_NEURALNETWORK_TRVAR  Evaluate training variance vector
%
%  Description
%    C = GP_NEURALNETWORK_TRVAR(GPCF, TX) takes in covariance
%    function of a Gaussian process GPCF and matrix TX that
%    contains training inputs. Returns variance vector C. Every
%    element i of C contains variance of input i in TX. This is 
%    a mandatory subfunction used for example in prediction and 
%    energy computations.
%
%
%  See also
%    GPCF_NEURALNETWORK_COV, GP_COV, GP_TRCOV
  
  if isfield(gpcf, 'selectedVariables')
    x=x(:,gpcf.selectedVariables); 
  end
  
  [n,m]=size(x);
  x_aug=[ones(n,1) x];
  
  if length(gpcf.weightSigma2) == 1
    % In the case of an isotropic NEURALNETWORK
    s = gpcf.weightSigma2*ones(1,m);
  else
    s = gpcf.weightSigma2;
  end

  s_tmp=sum(repmat([gpcf.biasSigma2 s], n, 1).*x_aug.^2,2);
  
  C=2/pi*asin(2*s_tmp./(1+2*s_tmp));
  C(C<eps)=0;
  
end

function reccf = gpcf_neuralnetwork_recappend(reccf, ri, gpcf)
%RECAPPEND Record append
%
%  Description
%    RECCF = GPCF_NEURALNETWORK_RECAPPEND(RECCF, RI, GPCF) takes
%    a covariance function record structure RECCF, record index
%    RI and covariance function structure GPCF with the current
%    MCMC samples of the parameters. Returns RECCF which contains
%    all the old samples and the current samples from GPCF. 
%    This subfunction is needed when using MCMC sampling (gp_mc).
%
%  See also
%    GP_MC and GP_MC -> RECAPPEND


  if nargin == 2
    % Initialize the record
    reccf.type = 'gpcf_neuralnetwork';
    reccf.nin = ri;
    reccf.nout = 1;

    % Initialize parameters
    reccf.weightSigma2= [];
    reccf.biasSigma2 = [];

    % Set the function handles
    reccf.fh.pak = @gpcf_neuralnetwork_pak;
    reccf.fh.unpak = @gpcf_neuralnetwork_unpak;
    reccf.fh.lp = @gpcf_neuralnetwork_lp;
    reccf.fh.lpg = @gpcf_neuralnetwork_lpg;
    reccf.fh.cfg = @gpcf_neuralnetwork_cfg;
    reccf.fh.cov = @gpcf_neuralnetwork_cov;
    reccf.fh.trcov  = @gpcf_neuralnetwork_trcov;
    reccf.fh.trvar  = @gpcf_neuralnetwork_trvar;
    reccf.fh.recappend = @gpcf_neuralnetwork_recappend;
    reccf.p=[];
    reccf.p.weightSigma2=[];
    reccf.p.biasSigma2=[];
    if ~isempty(ri.p.weightSigma2)
      reccf.p.weightSigma2 = ri.p.weightSigma2;
    end
    if ~isempty(ri.p.biasSigma2)
      reccf.p.biasSigma2 = ri.p.biasSigma2;
    end
  else
    % Append to the record
    gpp = gpcf.p;
    
    % record weightSigma2
    reccf.weightSigma2(ri,:)=gpcf.weightSigma2;
    if isfield(gpp,'weightSigma2') && ~isempty(gpp.weightSigma2)
      reccf.p.weightSigma2 = gpp.weightSigma2.fh.recappend(reccf.p.weightSigma2, ri, gpcf.p.weightSigma2);
    end
    
    % record biasSigma2
    reccf.biasSigma2(ri,:)=gpcf.biasSigma2;
    if isfield(gpp,'biasSigma2') && ~isempty(gpp.biasSigma2)
      reccf.p.biasSigma2 = gpp.biasSigma2.fh.recappend(reccf.p.biasSigma2, ri, gpcf.p.biasSigma2);
    end
  
    if isfield(gpcf, 'selectedVariables')
      reccf.selectedVariables = gpcf.selectedVariables;
    end
  end
end
