function gpcf = gpcf_rq(varargin)
%GPCF_RQ  Create a rational quadratic covariance function
%
%  Description
%    GPCF = GPCF_RQ('PARAM1',VALUE1,'PARAM2,VALUE2,...) creates
%    rational quadratic covariance function structure in which the
%    named parameters have the specified values. Any unspecified
%    parameters are set to default values.
%
%    GPCF = GPCF_RQ(GPCF,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a covariance function structure with the named
%    parameters altered with the specified values.
%  
%    Parameters for rational quadratic covariance function [default]
%      magnSigma2        - magnitude (squared) [0.1]
%      lengthScale       - length scale for each input. [1]
%                          This can be either scalar corresponding
%                          to an isotropic function or vector
%                          defining own length-scale for each input
%                          direction.
%      alpha             - shape parameter [20] 
%      magnSigma2_prior  - prior for magnSigma2  [prior_logunif]
%      lengthScale_prior - prior for lengthScale [prior_t]
%      alpha_prior       - prior for alpha [prior_unif]
%      metric            - metric structure used by the covariance function []
%      selectedVariables - vector defining which inputs are used [all]
%                          selectedVariables is shorthand for using
%                          metric_euclidean with corresponding components
%
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc.
%
%  See also
%    GP_SET, GPCF_*, PRIOR_*, METRIC_*

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2010 Tuomas Nikoskinen, Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  if nargin>0 && ischar(varargin{1}) && ismember(varargin{1},{'init' 'set'})
    % remove init and set
    varargin(1)=[];
  end
  
  ip=inputParser;
  ip.FunctionName = 'GPCF_RQ';
  ip.addOptional('gpcf', [], @isstruct);
  ip.addParamValue('magnSigma2',0.1, @(x) isscalar(x) && x>0);
  ip.addParamValue('lengthScale',1, @(x) isvector(x) && all(x>0));
  ip.addParamValue('alpha',20, @(x) isscalar(x) && x>0);
  ip.addParamValue('metric',[], @isstruct);
  ip.addParamValue('magnSigma2_prior', prior_logunif(), ...
                   @(x) isstruct(x) || isempty(x));
  ip.addParamValue('lengthScale_prior',prior_t(), ...
                   @(x) isstruct(x) || isempty(x));
  ip.addParamValue('alpha_prior', prior_unif(), ...
                   @(x) isstruct(x) || isempty(x));
  ip.addParamValue('selectedVariables',[], @(x) isempty(x) || ...
                   (isvector(x) && all(x>0)));
  ip.parse(varargin{:});
  gpcf=ip.Results.gpcf;

  if isempty(gpcf)
    init=true;
    gpcf.type = 'gpcf_rq';
  else
    if ~isfield(gpcf,'type') && ~isequal(gpcf.type,'gpcf_rq')
      error('First argument does not seem to be a valid covariance function structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('lengthScale',ip.UsingDefaults)
    gpcf.lengthScale = ip.Results.lengthScale;
  end
  if init || ~ismember('magnSigma2',ip.UsingDefaults)
    gpcf.magnSigma2 = ip.Results.magnSigma2;
  end
  if init || ~ismember('alpha',ip.UsingDefaults)
    gpcf.alpha = ip.Results.alpha;
  end

  % Initialize prior structure
  if init
    gpcf.p=[];
  end
  if init || ~ismember('lengthScale_prior',ip.UsingDefaults)
    gpcf.p.lengthScale=ip.Results.lengthScale_prior;
  end
  if init || ~ismember('magnSigma2_prior',ip.UsingDefaults)
    gpcf.p.magnSigma2=ip.Results.magnSigma2_prior;
  end
  if init || ~ismember('alpha_prior',ip.UsingDefaults)
    gpcf.p.alpha=ip.Results.alpha_prior;
  end

  %Initialize metric
  if ~ismember('metric',ip.UsingDefaults)
    if ~isempty(ip.Results.metric)
      gpcf.metric = ip.Results.metric;
      gpcf = rmfield(gpcf, 'lengthScale');
      gpcf.p = rmfield(gpcf.p, 'lengthScale');
    elseif isfield(gpcf,'metric')
      if ~isfield(gpcf,'lengthScale')
        gpcf.lengthScale = gpcf.metric.lengthScale;
      end
      if ~isfield(gpcf.p,'lengthScale')
        gpcf.p.lengthScale = gpcf.metric.p.lengthScale;
      end
      gpcf = rmfield(gpcf, 'metric');
    end
  end
  
  % selectedVariables options implemented using metric_euclidean
  if ~ismember('selectedVariables',ip.UsingDefaults)
    if ~isfield(gpcf,'metric')
      if ~isempty(ip.Results.selectedVariables)
        gpcf.selectedVariables = ip.Results.selectedVariables;
%         gpcf.metric=metric_euclidean('components',...
%                                      num2cell(ip.Results.selectedVariables),...
%                                      'lengthScale',gpcf.lengthScale,...
%                                      'lengthScale_prior',gpcf.p.lengthScale);
%         gpcf = rmfield(gpcf, 'lengthScale');
%         gpcf.p = rmfield(gpcf.p, 'lengthScale');
      end
    elseif isfield(gpcf,'metric') 
      if ~isempty(ip.Results.selectedVariables)
        gpcf.metric=metric_euclidean(gpcf.metric,...
                                     'components',...
                                     num2cell(ip.Results.selectedVariables));
        if ~ismember('lengthScale',ip.UsingDefaults)
          gpcf.metric.lengthScale=ip.Results.lengthScale;
          gpcf = rmfield(gpcf, 'lengthScale');
        end
        if ~ismember('lengthScale_prior',ip.UsingDefaults)
          gpcf.metric.p.lengthScale=ip.Results.lengthScale_prior;
          gpcf.p = rmfield(gpcf.p, 'lengthScale');
        end
      else
        if ~isfield(gpcf,'lengthScale')
          gpcf.lengthScale = gpcf.metric.lengthScale;
        end
        if ~isfield(gpcf.p,'lengthScale')
          gpcf.p.lengthScale = gpcf.metric.p.lengthScale;
        end
        gpcf = rmfield(gpcf, 'metric');
      end
    end
  end

  if init
    % Set the function handles to the subfunctions
    gpcf.fh.pak = @gpcf_rq_pak;
    gpcf.fh.unpak = @gpcf_rq_unpak;
    gpcf.fh.lp = @gpcf_rq_lp;
    gpcf.fh.lpg = @gpcf_rq_lpg;
    gpcf.fh.cfg = @gpcf_rq_cfg;
    gpcf.fh.ginput = @gpcf_rq_ginput;
    gpcf.fh.cov = @gpcf_rq_cov;
    gpcf.fh.trcov  = @gpcf_rq_trcov;
    gpcf.fh.trvar  = @gpcf_rq_trvar;
    gpcf.fh.recappend = @gpcf_rq_recappend;
  end

end

function [w, s] = gpcf_rq_pak(gpcf)
%GPCF_RQ_PAK  Combine GP covariance function parameters into
%             one vector
%
%  Description
%    W = GPCF_RQ_PAK(GPCF) takes a covariance function structure
%    GPCF and combines the covariance function parameters and
%    their hyperparameters into a single row vector W. This is a 
%    mandatory subfunction used for example in energy and gradient 
%    computations.
%
%       w = [ log(gpcf.magnSigma2)
%             (hyperparameters of gpcf.magnSigma2)
%             log(gpcf.lengthScale(:))
%             (hyperparameters of gpcf.lengthScale)
%             log(log(gpcf.alpha))
%             (hyperparameters of gpcf.alpha)]'
%
%  See also
%    GPCF_RQ_UNPAK

  w = []; s = {};
  
  if ~isempty(gpcf.p.magnSigma2)
    w = [w log(gpcf.magnSigma2)];
    s = [s; 'log(rq.magnSigma2)'];
    % Hyperparameters of magnSigma2
    [wh sh] = gpcf.p.magnSigma2.fh.pak(gpcf.p.magnSigma2);
    w = [w wh];
    s = [s; sh];
  end        

  if isfield(gpcf,'metric')
    [wm sm] = gpcf.metric.fh.pak(gpcf.metric);
    w = [w wm];
    s = [s; sm];
  else
    if ~isempty(gpcf.p.lengthScale)
      w = [w log(gpcf.lengthScale)];
      if numel(gpcf.lengthScale)>1
        s = [s; sprintf('log(rq.lengthScale x %d)',numel(gpcf.lengthScale))];
      else
        s = [s; 'log(rq.lengthScale)'];
      end
      % Hyperparameters of lengthScale
      [wh sh] = gpcf.p.lengthScale.fh.pak(gpcf.p.lengthScale);
      w = [w wh];
      s = [s; sh];
    end
  end
  
  if ~isempty(gpcf.p.alpha)
    w= [w log(log(gpcf.alpha))];
    s = [s; 'log(log(rq.alpha))'];
    % Hyperparameters of alpha
    [wh sh] = gpcf.p.alpha.fh.pak(gpcf.p.alpha);
    w = [w wh];
    s = [s; sh];
  end

end

function [gpcf, w] = gpcf_rq_unpak(gpcf, w)
%GPCF_RQ_UNPAK  Sets the covariance function parameters into
%                 the structure
%
%  Description
%    [GPCF, W] = GPCF_RQ_UNPAK(GPCF, W) takes a covariance
%    function structure GPCF and a hyper-parameter vector W, and
%    returns a covariance function structure identical to the
%    input, except that the covariance hyper-parameters have been
%    set to the values in W. Deletes the values set to GPCF from
%    W and returns the modified W. This is a mandatory subfunction 
%    used for example in energy and gradient computations.
%
%    Assignment is inverse of  
%       w = [ log(gpcf.magnSigma2)
%             (hyperparameters of gpcf.magnSigma2)
%             log(gpcf.lengthScale(:))
%             (hyperparameters of gpcf.lengthScale)
%             log(log(gpcf.alpha))
%             (hyperparameters of gpcf.alpha)]'
%
%  See also
%    GPCF_RQ_PAK
  
  gpp=gpcf.p;
  if ~isempty(gpp.magnSigma2)
    gpcf.magnSigma2 = exp(w(1));
    w = w(2:end);
    % Hyperparameters of magnSigma2
    [p, w] = gpcf.p.magnSigma2.fh.unpak(gpcf.p.magnSigma2, w);
    gpcf.p.magnSigma2 = p;
  end
  
  if isfield(gpcf,'metric')
    [metric, w] = gpcf.metric.fh.unpak(gpcf.metric, w);
    gpcf.metric = metric;
  else            
    if ~isempty(gpp.lengthScale)
      i1=1;
      i2=length(gpcf.lengthScale);
      gpcf.lengthScale = exp(w(i1:i2));
      w = w(i2+1:end);
      % Hyperparameters of lengthScale
      [p, w] = gpcf.p.lengthScale.fh.unpak(gpcf.p.lengthScale, w);
      gpcf.p.lengthScale = p;
    end
  end
  
  if ~isempty(gpp.alpha)
    gpcf.alpha = exp(exp(w(1)));
    w = w(2:end);
    % Hyperparameters of alpha
    [p, w] = gpcf.p.alpha.fh.unpak(gpcf.p.alpha, w);
    gpcf.p.alpha = p;
  end
  
end

function lp =gpcf_rq_lp(gpcf, x, t)
%GPCF_RQ_LP  Evaluate the log prior of covariance function parameters
%
%  Description
%    LP = GPCF_RQ_LP(GPCF, X, T) takes a covariance function
%    structure GPCF and returns log(p(th)), where th collects the
%    parameters. This is a mandatory subfunction used for example 
%    in energy computations.
%
%  See also
%    GPCF_RQ_PAK, GPCF_RQ_UNPAK, GPCF_RQ_LPG, GP_E

% Evaluate the prior contribution to the error. The parameters that
% are sampled are transformed, e.g., W = log(w) where w is all
% the "real" samples. On the other hand errors are evaluated in
% the W-space so we need take into account also the Jacobian of
% transformation, e.g., W -> w = exp(W). See Gelman et.al., 2004,
% Bayesian data Analysis, second edition, p24.
  lp = 0;
  gpp=gpcf.p;
  
  if ~isempty(gpcf.p.magnSigma2)
    lp = lp +gpp.magnSigma2.fh.lp(gpcf.magnSigma2, ...
                   gpp.magnSigma2) +log(gpcf.magnSigma2);
  end

  if isfield(gpcf,'metric')
    lp = lp +gpcf.metric.fh.lp(gpcf.metric);
  elseif ~isempty(gpp.lengthScale)
    lp = lp +gpp.lengthScale.fh.lp(gpcf.lengthScale, ...
                   gpp.lengthScale) +sum(log(gpcf.lengthScale));
  end

  if ~isempty(gpcf.p.alpha)
    lp = lp +gpp.alpha.fh.lp(gpcf.alpha, gpp.alpha) ...
         +log(gpcf.alpha) +log(log(gpcf.alpha));
  end
end

function lpg = gpcf_rq_lpg(gpcf)
%GPCF_RQ_LPG  Evaluate gradient of the log prior with respect
%             to the parameters
%
%  Description
%    LPG = GPCF_RQ_LPG(GPCF) takes a covariance function
%    structure GPCF and returns LPG = d log (p(th))/dth, where th
%    is the vector of parameters. This is a mandatory subfunction 
%    used for example in energy and gradient computations.
%
%  See also
%    GPCF_RQ_PAK, GPCF_RQ_UNPAK, GPCF_RQ_LP, GP_G

  lpg = [];
  gpp=gpcf.p;
  
  if ~isempty(gpcf.p.magnSigma2)            
    lpgs = gpp.magnSigma2.fh.lpg(gpcf.magnSigma2, gpp.magnSigma2);
    lpg = [lpg lpgs(1).*gpcf.magnSigma2+1 lpgs(2:end)];
  end
  
  if isfield(gpcf,'metric')
    lpg_dist = gpcf.metric.fh.lpg(gpcf.metric);
    lpg=[lpg lpg_dist];
  else
    if ~isempty(gpcf.p.lengthScale)
      lll = length(gpcf.lengthScale);
      lpgs = gpp.lengthScale.fh.lpg(gpcf.lengthScale, gpp.lengthScale);
      lpg = [lpg lpgs(1:lll).*gpcf.lengthScale+1 lpgs(lll+1:end)];
    end
  end
  
  if ~isempty(gpcf.p.alpha)            
    lpgs = gpp.alpha.fh.lpg(gpcf.alpha, gpp.alpha);
    lpg = [lpg lpgs(1).*gpcf.alpha.*log(gpcf.alpha)+log(gpcf.alpha)+1 lpgs(2:end)];
  end
  
end

function DKff = gpcf_rq_cfg(gpcf, x, x2, mask, i1)
%GPCF_RQ_CFG  Evaluate gradient of covariance function
%             with respect to the parameters
%
%  Description
%    DKff = GPCF_RQ_CFG(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to th (cell array with matrix elements). This is a 
%    mandatory subfunction used in gradient computations.
%
%    DKff = GPCF_RQ_CFG(GPCF, X, X2) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X2) with
%    respect to th (cell array with matrix elements). This subfunction 
%    is needed when using sparse approximations (e.g. FIC).
%
%    DKff = GPCF_RQ_CFG(GPCF, X, [], MASK) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the diagonal of gradients of covariance matrix
%    Kff = k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse
%    approximations (e.g. FIC).
%
%    DKff = GPCF_RQ_CFG(GPCF, X, X2, [], i) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X2), or 
%    k(X,X) if X2 is empty, with respect to ith hyperparameter. This
%    subfunction is needed when using memory save option in gp_set.
%
%  See also
%   GPCF_RQ_PAK, GPCF_RQ_UNPAK, GPCF_RQ_LP, GP_G

  gpp=gpcf.p;
  a=(gpcf.alpha+1)/gpcf.alpha;

  i2=1;
  DKff = {};

  if nargin==5
    % Use memory save option
    savememory=1;
    if i1==0
      % Return number of hyperparameters
      i=0;
      if ~isempty(gpcf.p.magnSigma2)
        i=i+1;
      end
      if ~isempty(gpcf.p.lengthScale)
        i=i+length(gpcf.lengthScale);
      end
      if ~isempty(gpcf.p.alpha)
        i=i+1;
      end
      DKff=i;
      return
    end
  else
    savememory=0;
  end

  % Evaluate: DKff{1} = d Kff / d magnSigma2
  %           DKff{2} = d Kff / d alpha
  %           DKff{3} = d Kff / d lengthscale
  % NOTE! Here we have already taken into account that the parameters
  % are transformed through log() and thus dK/dlog(p) = p * dK/dp
  % (or loglog gor alpha)

  % evaluate the gradient for training covariance
  if nargin == 2 || (isempty(x2) && isempty(mask))
    Cdm = gpcf_rq_trcov(gpcf, x);
    ii1=0;

    if ~isempty(gpcf.p.magnSigma2)
      ii1 = ii1 +1;
      DKff{ii1} = Cdm;
    end
    
    ma2=gpcf.magnSigma2;
    
    if isfield(gpcf,'metric')
      dist = gpcf.metric.fh.dist(gpcf.metric, x);
      distg = gpcf.metric.fh.distg(gpcf.metric, x);
      gprior_dist = gpcf.metric.fh.lpg(gpcf.metric);
      % dalpha
      ii1=ii1+1;
      DKff{ii1} = (ma2.^(1-a).*.5.*dist.^2.*Cdm.^a - gpcf.alpha.*log(Cdm.^(-1/gpcf.alpha)./ma2.^(-1/gpcf.alpha)).*Cdm).*log(gpcf.alpha);

      % dlengthscale
      for i=1:length(distg)
        ii1=ii1+1;
        DKff{ii1} = Cdm.*-dist./(1+dist.^2./(2*gpcf.alpha)).*distg{i};
      end
    else
      if isfield(gpcf, 'selectedVariables')
        x = x(:,gpcf.selectedVariables);
      end
      [n, m] =size(x);
      if ~savememory
        i1=1:m;
        i2=1:m;
      else
        if i1==1
          DKff=DKff{1};
          return
        end
        i1=i1-1;
        if i1 > length(gpcf.lengthScale)
          i2=1:m;
        else
          i2=i1;
        end
        ii1=ii1-1;
      end
      % loop over all the lengthScales
      if length(gpcf.lengthScale) == 1
        % Isotropic = no ARD
        s = 1./(gpcf.lengthScale^2);
        dist2 = 0;
        for i=1:m
          dist2 = dist2 + (bsxfun(@minus,x(:,i),x(:,i)')).^2;
        end
        if ~isempty(gpcf.p.lengthScale) && (~savememory || i1==1)
          % dlengthscale
          ii1 = ii1+1;
          DKff{ii1} = Cdm.^a.*s.*dist2.*gpcf.magnSigma2^(-a+1);
        end
        if ~isempty(gpcf.p.alpha) && (~savememory || length(DKff) == 1)
          % dalpha
          ii1=ii1+1;
          DKff{ii1} = (ma2^(1-a).*.5.*dist2.*s.*Cdm.^a - gpcf.alpha.*log(Cdm.^(-1/gpcf.alpha)./ma2^(-1/gpcf.alpha)).*Cdm).*log(gpcf.alpha);
        end
      else
        % ARD
        s = 1./(gpcf.lengthScale.^2);
        D=zeros(size(Cdm));
        for i=i2
          dist2 =(bsxfun(@minus,x(:,i),x(:,i)')).^2;
          % sum distance for the dalpha
          D=D+dist2.*s(i); 
          % dlengthscale
          if ~isempty(gpcf.p.lengthScale) && all(i1 <= m)
            ii1 = ii1+1;
            DKff{ii1}=Cdm.^a.*s(i).*dist2.*gpcf.magnSigma2.^(-a+1);
          end
        end
        if ~isempty(gpcf.p.alpha) && (~savememory || isvector(i2))
          % dalpha
          ii1=ii1+1;
          DKff{ii1} = (ma2^(1-a).*.5.*D.*Cdm.^a - gpcf.alpha.*log(Cdm.^(-1/gpcf.alpha)./ma2^(-1/gpcf.alpha)).*Cdm).*log(gpcf.alpha);
        end
      end
    end
    % Evaluate the gradient of non-symmetric covariance (e.g. K_fu)
  elseif nargin == 3 || isempty(mask)
    if size(x,2) ~= size(x2,2)
      error('gpcf_rq -> _ghyper: The number of columns in x and x2 has to be the same. ')
    end
    
    ii1=0;
    K = gpcf.fh.cov(gpcf, x, x2);
    if ~isempty(gpcf.p.magnSigma2)
      ii1=ii1+1;
      DKff{ii1} = K;
    end
    
    if isfield(gpcf,'metric')                
      dist = gpcf.metric.fh.dist(gpcf.metric, x, x2);
      distg = gpcf.metric.fh.distg(gpcf.metric, x, x2);
      gprior_dist = gpcf.metric.fh.lpg(gpcf.metric);
      for i=1:length(distg)
        ii1 = ii1+1;                    
        DKff{ii1} = -K.*distg{i};                    
      end
    else
      if isfield(gpcf, 'selectedVariables')
        x = x(:,gpcf.selectedVariables);
        x2 = x2(:,gpcf.selectedVariables);
      end
      [n, m] =size(x);
      if ~savememory
        i1=1:m;
        i2=1:m;
      else
        if i1==1
          DKff=DKff{1};
          return
        end
        i1=i1-1;
        if i1 > length(gpcf.lengthScale)
          i2=1:m;
        else
          i2=i1;
        end
        ii1=ii1-1;
      end
      % Evaluate help matrix for calculations of derivatives with respect to the lengthScale
      if length(gpcf.lengthScale) == 1
        % In the case of an isotropic EXP
        s = 1/gpcf.lengthScale^2;
        dist = 0;
        for i=1:m
          dist = dist + (bsxfun(@minus,x(:,i),x2(:,i)')).^2;
        end
        DK_l = s.*K.^a.*dist.*gpcf.magnSigma2^(1-a);
        ii1=ii1+1;
        DKff{ii1} = DK_l;
        if ~isempty(gpcf.p.alpha) && (~savememory || length(DKff) == 1)
          % dalpha
          ii1=ii1+(1-savememory);
          DKff{ii1} = (gpcf.magnSigma2^(1-a).*.5.*dist.*s.*K.^a - gpcf.alpha.*log(K.^(-1/gpcf.alpha)./gpcf.magnSigma2^(-1/gpcf.alpha)).*K).*log(gpcf.alpha);
        end
      else
        % In the case ARD is used
        s = 1./gpcf.lengthScale.^2;        % set the length
        D=zeros(size(K));
        for i=i2
          dist2 =(bsxfun(@minus,x(:,i),x2(:,i)')).^2;
          % sum distance for the dalpha
          D=D+dist2.*s(i);
          if ~isempty(gpcf.p.lengthScale) && all(i1 <= m)
            D1 = s(i).*K.^a.*dist2.*gpcf.magnSigma2^(1-a);
            ii1=ii1+1;
            DKff{ii1} = D1;
          end
        end
        if ~isempty(gpcf.p.alpha) && (~savememory || isvector(i2))
          % dalpha
          ii1=ii1+1;
          DKff{ii1} = (gpcf.magnSigma2^(1-a).*.5.*D.*K.^a - gpcf.alpha.*log(K.^(-1/gpcf.alpha)./gpcf.magnSigma2^(-1/gpcf.alpha)).*K).*log(gpcf.alpha);
        end
      end
    end
    % Evaluate: DKff{1}    = d mask(Kff,I) / d magnSigma2
    %           DKff{2...} = d mask(Kff,I) / d lengthScale
  elseif nargin == 4 || nargin == 5
    if isfield(gpcf,'metric')
      ii1=1;
      [n, m] =size(x);
      DKff{ii1} = gpcf.fh.trvar(gpcf, x);   % d mask(Kff,I) / d magnSigma2
      
      dist = 0;
      distg = gpcf.metric.fh.distg(gpcf.metric, x, [], 1);
      gprior_dist = gpcf.metric.fh.lpg(gpcf.metric);
      for i=1:length(distg)
        ii1 = ii1+1;
        DKff{ii1} = 0;
      end
    else
      ii1=0;
      if ~isempty(gpcf.p.magnSigma2) && (~savememory || all(i1==1))
        ii1=ii1+1;
        DKff{ii1} = gpcf.fh.trvar(gpcf, x);   % d mask(Kff,I) / d magnSigma2
      end
      for i2=1:length(gpcf.lengthScale)
        ii1 = ii1+1;
        DKff{ii1}  = 0;                          % d mask(Kff,I) / d lengthScale
      end
    end
  end
  if savememory
    DKff=DKff{1};
  end
end

function DKff = gpcf_rq_ginput(gpcf, x, x2, i1)
%GPCF_RQ_GINPUT  Evaluate gradient of covariance function with 
%                respect to x
%
%  Description
%    DKff = GPCF_RQ_GINPUT(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to X (cell array with matrix elements). This 
%    subfunction is needed when computing gradients with respect 
%    to inducing inputs in sparse approximations.
%
%    DKff = GPCF_RQ_GINPUT(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors
%    and returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_RQ_GINPUT(GPCF, X, X2, i) takes a covariance
%    function structure GPCF, a matrix X of input vectors
%    and returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2), or k(X,X) if X2 is empty, with respect to ith 
%    covariate in X (cell array with matrix elements). This
%    subfunction is needed when using memory save option in 
%    gp_set.
%
%  See also
%   GPCF_RQ_PAK, GPCF_RQ_UNPAK, GPCF_RQ_LP, GP_G
  
  a=(gpcf.alpha+1)/gpcf.alpha;
  [n, m] =size(x);
  
  if nargin<4
    i1=1:m;
  else
    % Use memory save option
    if i1==0
      % Return number of covariates
      if isfield(gpcf,'selectedVariables')
        DKff=length(gpcf.selectedVariables);
      else
        DKff=m;
      end
      return
    end
  end
  if nargin == 2 || isempty(x2)
    K = gpcf.fh.trcov(gpcf, x);
    ii1 = 0;
    if isfield(gpcf,'metric')
      dist = gpcf.metric.fh.dist(gpcf.metric, x);
      [gdist, gprior_dist] = gpcf.metric.fh.ginput(gpcf.metric, x);
      for i=1:length(gdist)
        ii1 = ii1+1;
        DKff{ii1} = -K.*gdist{ii1};
        gprior(ii1) = gprior_dist(ii1);
      end
    else
      if length(gpcf.lengthScale) == 1
        % In the case of an isotropic RQ
        s = repmat(1./gpcf.lengthScale.^2, 1, m);
      else
        s = 1./gpcf.lengthScale.^2;
      end
      for i=i1
        for j = 1:n
          DK = zeros(size(K));
          DK(j,:) = -s(i).*bsxfun(@minus,x(j,i),x(:,i)');
          DK = DK + DK';    
          
          DK = DK.*K.^a.*gpcf.magnSigma2^(1-a);      
          
          ii1 = ii1 + 1;
          DKff{ii1} = DK;
          gprior(ii1) = 0; 
        end
      end
    end
    
  elseif nargin == 3 || nargin == 4
    [n2, m2] =size(x2);
    K = gpcf.fh.cov(gpcf, x, x2);
    ii1 = 0;
    if isfield(gpcf,'metric')
      dist = gpcf.metric.fh.dist(gpcf.metric, x, x2);
      [gdist, gprior_dist] = gpcf.metric.fh.ginput(gpcf.metric, x, x2);
      for i=1:length(gdist)
        ii1 = ii1+1;
        DKff{ii1}   = -K.*gdist{ii1};
        gprior(ii1) = gprior_dist(ii1);
      end
    else 
      if length(gpcf.lengthScale) == 1
        % In the case of an isotropic RQ
        s = repmat(1./gpcf.lengthScale.^2, 1, m);
      else
        s = 1./gpcf.lengthScale.^2;
      end
      
      ii1 = 0;
      for i=i1
        for j = 1:n
          DK= zeros(size(K));
          DK(j,:) = -s(i).*bsxfun(@minus,x(j,i),x2(:,i)');
          
          DK = DK.*K.^a.*gpcf.magnSigma2^(1-a);
          
          ii1 = ii1 + 1;
          DKff{ii1} = DK;
          gprior(ii1) = 0; 
        end
      end
    end
  end
end

function C = gpcf_rq_cov(gpcf, x1, x2)
% GP_RQ_COV  Evaluate covariance matrix between two input vectors
%
%  Description         
%    C = GP_RQ_COV(GP, TX, X) takes in covariance function of a
%    Gaussian process GP and two matrixes TX and X that contain
%    input vectors to GP. Returns covariance matrix C. Every
%    element ij of C contains covariance between inputs i in TX
%    and j in X. This is a mandatory subfunction used for example in
%    prediction and energy computations.
%
%  See also
%    GPCF_RQ_TRCOV, GPCF_RQ_TRVAR, GP_COV, GP_TRCOV
  
  if isempty(x2)
    x2=x1;
  end

  if size(x1,2)~=size(x2,2)
    error('the number of columns of X1 and X2 has to be same')
  end

  if isfield(gpcf,'metric')
    dist = gpcf.metric.fh.dist(gpcf.metric, x1, x2).^2;
    dist(dist<eps) = 0;
    C = gpcf.magnSigma2.*(1+dist./(2*gpcf.alpha)).^(-gpcf.alpha);
  else
    if isfield(gpcf, 'selectedVariables')
      x1 = x1(:,gpcf.selectedVariables);
      x2 = x2(:,gpcf.selectedVariables);
    end
    [n1,m1]=size(x1);
    [n2,m2]=size(x2);
    C=zeros(n1,n2);
    ma2 = gpcf.magnSigma2;
    
    % Evaluate the covariance
    if ~isempty(gpcf.lengthScale)  
      s2 = 1./(2.*gpcf.alpha.*gpcf.lengthScale.^2);      
      % If ARD is not used make s a vector of 
      % equal elements 
      if size(s2)==1
        s2 = repmat(s2,1,m1);
      end
      dist=zeros(n1,n2);
      for j=1:m1
        dist = dist + s2(j).*(bsxfun(@minus,x1(:,j),x2(:,j)')).^2;
      end
      dist(dist<eps) = 0;
      C = ma2.*(1+dist).^(-gpcf.alpha);
    end
  end
end

function C = gpcf_rq_trcov(gpcf, x)
%GP_RQ_TRCOV  Evaluate training covariance matrix of inputs
%
%  Description
%    C = GP_RQ_TRCOV(GP, TX) takes in covariance function of a
%    Gaussian process GP and matrix TX that contains training
%    input vectors. Returns covariance matrix C. Every element ij
%    of C contains covariance between inputs i and j in TX. This 
%    is a mandatory subfunction used for example in prediction and 
%    energy computations.
%
%  See also
%    GPCF_RQ_COV, GPCF_RQ_TRVAR, GP_COV, GP_TRCOV

  if isfield(gpcf,'metric')
    % If other than scaled euclidean metric
    [n, m] =size(x);            
    ma2 = gpcf.magnSigma2;
    
    C = zeros(n,n);
    for ii1=1:n-1
      d = zeros(n-ii1,1);
      col_ind = ii1+1:n;
      d = gpcf.metric.fh.dist(gpcf.metric, x(col_ind,:), x(ii1,:)).^2;                
      C(col_ind,ii1) = d;
    end
    C(C<eps) = 0;
    C = C+C';
    C = ma2.*(1+C./(2*gpcf.alpha)).^(-gpcf.alpha);     
  else
    % If scaled euclidean metric
    % Try to use the C-implementation
    C=trcov(gpcf, x);
    if isnan(C)
      % If there wasn't C-implementation do here
      if isfield(gpcf, 'selectedVariables')
        x = x(:,gpcf.selectedVariables);
      end
      [n, m] =size(x);
      
      s2 = 1./(2*gpcf.alpha.*gpcf.lengthScale.^2);
      if size(s2)==1
        s2 = repmat(s2,1,m);
      end
      ma2 = gpcf.magnSigma2;
      
      C = zeros(n,n);
      for ii1=1:n-1
        d = zeros(n-ii1,1);
        col_ind = ii1+1:n;
        for ii2=1:m
          d = d+s2(ii2).*(x(col_ind,ii2)-x(ii1,ii2)).^2;
        end
        C(col_ind,ii1) = d;
      end
      C(C<eps) = 0;
      C = C+C';
      C = ma2.*(1+C).^(-gpcf.alpha);
    end
  end
end

function C = gpcf_rq_trvar(gpcf, x)
%GP_RQ_TRVAR  Evaluate training variance vector
%
%  Description
%    C = GP_RQ_TRVAR(GPCF, TX) takes in covariance function of a
%    Gaussian process GPCF and matrix TX that contains training
%    inputs. Returns variance vector C. Every element i of C
%    contains variance of input i in TX. This is a mandatory 
%    subfunction used for example in prediction and energy computations.
%
%  See also
%    GPCF_RQ_COV, GP_COV, GP_TRCOV

  [n, m] =size(x);

  C = ones(n,1).*gpcf.magnSigma2;
  C(C<eps)=0;
end

function reccf = gpcf_rq_recappend(reccf, ri, gpcf)
%RECAPPEND  Record append
%
%  Description
%    RECCF = GPCF_RQ_RECAPPEND(RECCF, RI, GPCF) takes a
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
    reccf.type = 'gpcf_rq';

    % Initialize parameters
    reccf.lengthScale= [];
    reccf.magnSigma2 = [];
    reccf.gpcf.alpha = [];
    
    % Set the function handles
    reccf.fh.pak = @gpcf_rq_pak;
    reccf.fh.unpak = @gpcf_rq_unpak;
    reccf.fh.e = @gpcf_rq_lp;
    reccf.fh.g = @gpcf_rq_g;
    reccf.fh.cov = @gpcf_rq_cov;
    reccf.fh.trcov  = @gpcf_rq_trcov;
    reccf.fh.trvar  = @gpcf_rq_trvar;
    reccf.fh.recappend = @gpcf_rq_recappend;  
    reccf.p=[];
    reccf.p.lengthScale=[];
    reccf.p.magnSigma2=[];
    if isfield(ri.p,'lengthScale') && ~isempty(ri.p.lengthScale)
      reccf.p.lengthScale = ri.p.lengthScale;
    end
    if ~isempty(ri.p.magnSigma2)
      reccf.p.magnSigma2 = ri.p.magnSigma2;
    end
    if ~isempty(ri.p.alpha)
      reccf.p.alpha = ri.p.alpha;
    end
  else
    % Append to the record

  gpp = gpcf.p;

    if ~isfield(gpcf,'metric')
      % record lengthScale
      reccf.lengthScale(ri,:)=gpcf.lengthScale;
      if isfield(gpp,'lengthScale') && ~isempty(gpp.lengthScale)
        reccf.p.lengthScale = gpp.lengthScale.fh.recappend(reccf.p.lengthScale, ri, gpcf.p.lengthScale);
      end
    end
    
    % record magnSigma2
    reccf.magnSigma2(ri,:)=gpcf.magnSigma2;
    if isfield(gpp,'magnSigma2') && ~isempty(gpp.magnSigma2)
      reccf.p.magnSigma2 = gpp.magnSigma2.fh.recappend(reccf.p.magnSigma2, ri, gpcf.p.magnSigma2);
    end

    reccf.alpha(ri,:)=gpcf.alpha;
    if isfield(gpp,'alpha') && ~isempty(ri.p.alpha)
      reccf.p.alpha = gpp.alpha.fh.recappend(reccf.p.alpha, ri, gpcf.p.alpha);
    end
    
  end
end
