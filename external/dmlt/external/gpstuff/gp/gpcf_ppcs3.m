function gpcf = gpcf_ppcs3(varargin)
%GPCF_PPCS3  Create a piece wise polynomial (q=3) covariance function 
%
%  Description
%    GPCF = GPCF_PPCS3('nin',nin,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates piece wise polynomial (q=3) covariance function
%    structure in which the named parameters have the specified
%    values. Any unspecified parameters are set to default values. 
%    Obligatory parameter is 'nin', which tells the dimension
%    of input space.
%  
%    GPCF = GPCF_PPCS3(GPCF,'PARAM1',VALUE1,'PARAM2,VALUE2,...)
%    modify a covariance function structure with the named
%    parameters altered with the specified values.
%
%    Parameters for piece wise polynomial (q=3) covariance function [default]
%      magnSigma2        - magnitude (squared) [0.1]
%      lengthScale       - length scale for each input. [1]
%                          This can be either scalar corresponding
%                          to an isotropic function or vector
%                          defining own length-scale for each input
%                          direction.
%      l_nin             - order of the polynomial [floor(nin/2) + 4]
%                          Has to be greater than or equal to default.
%      magnSigma2_prior  - prior for magnSigma2  [prior_logunif]
%      lengthScale_prior - prior for lengthScale [prior_t]
%      metric            - metric structure used by the covariance function []
%      selectedVariables - vector defining which inputs are used [all]
%                          selectedVariables is shorthand for using
%                          metric_euclidean with corresponding components
%
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc.
%
%    The piecewise polynomial function is the following:
%
%      k(x_i, x_j) = ma2*cs^(l+3)*((l^3 + 9*l^2 + 23*l + 15)*r^3 + 
%                    (6*l^2 + 36*l + 45)*r^2 + (15*l + 45)*r + 15)/15
%
%      where r = sum( (x_i,d - x_j,d)^2/l^2_d )
%            l = floor(l_nin/2) + 4     
%            cs = max(0,1-r)
%      and l_nin must be greater or equal to gpcf.nin
%       
%    NOTE! Use of gpcf_ppcs3 requires that you have installed
%    GPstuff with SuiteSparse.
%
%  See also
%    GP_SET, GPCF_*, PRIOR_*, METRIC_*

% Copyright (c) 2009-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  if nargin>0 && ischar(varargin{1}) && ismember(varargin{1},{'init' 'set'})
    % remove init and set
    varargin(1)=[];
  end
  
  ip=inputParser;
  ip.FunctionName = 'GPCF_PPCS3';
  ip.addOptional('gpcf', [], @isstruct);
  ip.addParamValue('nin',[], @(x) isscalar(x) && x>0 && mod(x,1)==0);
  ip.addParamValue('magnSigma2',0.1, @(x) isscalar(x) && x>0);
  ip.addParamValue('lengthScale',1, @(x) isvector(x) && all(x>0));
  ip.addParamValue('l_nin',[], @(x) isscalar(x) && x>0 && mod(x,1)==0);
  ip.addParamValue('metric',[], @isstruct);
  ip.addParamValue('magnSigma2_prior', prior_logunif(), ...
                   @(x) isstruct(x) || isempty(x));
  ip.addParamValue('lengthScale_prior',prior_t(), ...
                   @(x) isstruct(x) || isempty(x));
  ip.addParamValue('selectedVariables',[], @(x) isempty(x) || ...
                   (isvector(x) && all(x>0)));
  ip.parse(varargin{:});
  gpcf=ip.Results.gpcf;

  if isempty(gpcf)
    % Check that SuiteSparse is available
    if ~exist('ldlchol')
      error('SuiteSparse is not installed (or it is not in the path). gpcf_ppcs3 cannot be used!')
    end
    init=true;
    gpcf.nin=ip.Results.nin;
    if isempty(gpcf.nin)
      error('nin has to be given for ppcs: gpcf_ppcs3(''nin'',NIN,...)')
    end
    gpcf.type = 'gpcf_ppcs3';
    % cf is compactly supported
    gpcf.cs = 1;
  else
    if ~isfield(gpcf,'type') && ~isequal(gpcf.type,'gpcf_ppcs3')
      error('First argument does not seem to be a valid covariance function structure')
    end
    init=false;
  end
  if init
    % Set the function handles to the subfunctions
    gpcf.fh.pak = @gpcf_ppcs3_pak;
    gpcf.fh.unpak = @gpcf_ppcs3_unpak;
    gpcf.fh.lp = @gpcf_ppcs3_lp;
    gpcf.fh.lpg = @gpcf_ppcs3_lpg;
    gpcf.fh.cfg = @gpcf_ppcs3_cfg;
    gpcf.fh.ginput = @gpcf_ppcs3_ginput;
    gpcf.fh.cov = @gpcf_ppcs3_cov;
    gpcf.fh.trcov  = @gpcf_ppcs3_trcov;
    gpcf.fh.trvar  = @gpcf_ppcs3_trvar;
    gpcf.fh.recappend = @gpcf_ppcs3_recappend;
  end

  % Initialize parameters
  if init || ~ismember('l_nin',ip.UsingDefaults)
    gpcf.l=ip.Results.l_nin;
    if isempty(gpcf.l)
      gpcf.l = floor(gpcf.nin/2) + 4;
    end
    if gpcf.l < gpcf.nin
      error('The l_nin has to be greater than or equal to the number of inputs!')
    end
  end
  if init || ~ismember('lengthScale',ip.UsingDefaults)
    gpcf.lengthScale = ip.Results.lengthScale;
  end
  if init || ~ismember('magnSigma2',ip.UsingDefaults)
    gpcf.magnSigma2 = ip.Results.magnSigma2;
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

end

function [w,s] = gpcf_ppcs3_pak(gpcf)
%GPCF_PPCS3_PAK  Combine GP covariance function parameters into
%                one vector
%
%  Description
%    W = GPCF_PPCS3_PAK(GPCF) takes a covariance function
%    structure GPCF and combines the covariance function
%    parameters and their hyperparameters into a single row
%    vector W. This is a mandatory subfunction used for 
%    example in energy and gradient computations.
%
%       w = [ log(gpcf.magnSigma2)
%             (hyperparameters of gpcf.magnSigma2)
%             log(gpcf.lengthScale(:))
%             (hyperparameters of gpcf.lengthScale)]'
%
%  See also
%    GPCF_PPCS3_UNPAK

  w = []; s = {};
  
  if ~isempty(gpcf.p.magnSigma2)
    w = [w log(gpcf.magnSigma2)];
    s = [s; 'log(ppcs3.magnSigma2)'];
    % Hyperparameters of magnSigma2
    [wh sh] = gpcf.p.magnSigma2.fh.pak(gpcf.p.magnSigma2);
    w = [w wh];
    s = [s; sh];
  end        

  if isfield(gpcf,'metric')
    [wh sh]=gpcf.metric.fh.pak(gpcf.metric);
    w = [w wh];
    s = [s; sh];
  else
    if ~isempty(gpcf.p.lengthScale)
      w = [w log(gpcf.lengthScale)];
      if numel(gpcf.lengthScale)>1
        s = [s; sprintf('log(ppcs3.lengthScale x %d)',numel(gpcf.lengthScale))];
      else
        s = [s; 'log(ppcs3.lengthScale)'];
      end
      % Hyperparameters of lengthScale
      [wh  sh] = gpcf.p.lengthScale.fh.pak(gpcf.p.lengthScale);
      w = [w wh];
      s = [s; sh];
    end
  end

end

function [gpcf, w] = gpcf_ppcs3_unpak(gpcf, w)
%GPCF_PPCS3_UNPAK  Sets the covariance function parameters into
%                  the structure
%
%  Description
%    [GPCF, W] = GPCF_PPCS3_UNPAK(GPCF, W) takes a covariance
%    function structure GPCF and a hyper-parameter vector W,
%    and returns a covariance function structure identical
%    to the input, except that the covariance hyper-parameters
%    have been set to the values in W. Deletes the values set to
%    GPCF from W and returns the modified W. This is a mandatory 
%    subfunction used for example in energy and gradient computations.
%
%    Assignment is inverse of  
%       w = [ log(gpcf.magnSigma2)
%             (hyperparameters of gpcf.magnSigma2)
%             log(gpcf.lengthScale(:))
%             (hyperparameters of gpcf.lengthScale)]'
%
%  See also
%    GPCF_PPCS3_PAK

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
  
end

function lp = gpcf_ppcs3_lp(gpcf)
%GPCF_PPCS3_LP  Evaluate the log prior of covariance function parameters
%
%  Description
%    LP = GPCF_PPCS3_LP(GPCF, X, T) takes a covariance function
%    structure GPCF and returns log(p(th)), where th collects the
%    parameters. This is a mandatory subfunction used for example 
%    in energy computations.
%
%  See also
%    GPCF_PPCS3_PAK, GPCF_PPCS3_UNPAK, GPCF_PPCS3_LPG, GP_E

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
end

function lpg = gpcf_ppcs3_lpg(gpcf)
%GPCF_PPCS3_LPG  Evaluate gradient of the log prior with respect
%               to the parameters.
%
%  Description
%    LPG = GPCF_PPCS3_LPG(GPCF) takes a covariance function
%    structure GPCF and returns LPG = d log (p(th))/dth, where th
%    is the vector of parameters. This is a mandatory subfunction 
%    used for example in gradient computations.
%
%  See also
%    GPCF_PPCS3_PAK, GPCF_PPCS3_UNPAK, GPCF_PPCS3_LP, GP_G

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
end

function DKff = gpcf_ppcs3_cfg(gpcf, x, x2, mask, i1)
%GPCF_PPCS3_CFG  Evaluate gradient of covariance function
%                with respect to the parameters
%
%  Description
%    DKff = GPCF_PPCS3_CFG(GPCF, X) takes a covariance function
%    structure GPCF, a matrix X of input vectors and returns
%    DKff, the gradients of covariance matrix Kff = k(X,X) with
%    respect to th (cell array with matrix elements). This is a 
%    mandatory subfunction used in gradient computations.
%
%    DKff = GPCF_PPCS3_CFG(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
%
%    DKff = GPCF_PPCS3_CFG(GPCF, X, X2, [], i) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2), or k(X,X) if X2 is empty, with respect to ith 
%    hyperparameter. This subfunction is needed when using sparse 
%    approximations (e.g. FIC).
% 
%    DKff = GPCF_PPCS3_CFG(GPCF, X, [], MASK) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the diagonal of gradients of covariance matrix
%    Kff = k(X,X2) with respect to th (cell array with matrix
%    elements). This subfunction is needed when using memory
%    save option in gp_set.
%
%  See also
%   GPCF_PPCS3_PAK, GPCF_PPCS3_UNPAK, GPCF_PPCS3_LP, GP_G

  gpp=gpcf.p;

  i2=1;
  DKff = {};
  gprior = [];
  
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
      DKff=i;
      return
    end
  else
    savememory=0;
  end

  % Evaluate: DKff{1} = d Kff / d magnSigma2
  %           DKff{2} = d Kff / d lengthScale
  % NOTE! Here we have already taken into account that the parameters
  % are transformed through log() and thus dK/dlog(p) = p * dK/dp

  % evaluate the gradient for training covariance
  if nargin == 2 || (isempty(x2) && isempty(mask))
    Cdm = gpcf_ppcs3_trcov(gpcf, x);

    ii1=0;

    if ~isempty(gpcf.p.magnSigma2)
      ii1 = ii1 +1;
      DKff{ii1} = Cdm;
    end
    
    l = gpcf.l;
    [I,J] = find(Cdm);
    
    if isfield(gpcf,'metric')
      % Compute the sparse distance matrix and its gradient.
      [n, m] =size(x);
      ntriplets = (nnz(Cdm)-n)./2;
      I = zeros(ntriplets,1);
      J = zeros(ntriplets,1);
      dist = zeros(ntriplets,1);
      for jj = 1:length(gpcf.metric.components)
        gdist{jj} = zeros(ntriplets,1);
      end
      ntriplets = 0;                
      for ii=1:n-1
        col_ind = ii + find(Cdm(ii+1:n,ii));
        d = gpcf.metric.fh.dist(gpcf.metric, x(col_ind,:), x(ii,:));
        gd = gpcf.metric.fh.distg(gpcf.metric, x(col_ind,:), x(ii,:));

        ntrip_prev = ntriplets;
        ntriplets = ntriplets + length(d);
        
        ind_tr = ntrip_prev+1:ntriplets;
        I(ind_tr) = col_ind;
        J(ind_tr) = ii;
        dist(ind_tr) = d;
        for jj = 1:length(gd)
          gdist{jj}(ind_tr) = gd{jj};
        end
      end
      
      ma2 = gpcf.magnSigma2;
      
      cs = 1-dist;
      
      const1 = l^3 + 9*l^2 + 23*l + 15;
      const2 = 6*l^2 + 36*l + 45;
      const3 = 15*l + 45;
      
      Dd = -(l+3).*cs.^(l+2).*(const1.*dist.^3 + const2.*dist.^2 + const3.*dist + 15)/15;
      Dd = Dd + cs.^(l+3).*(3.*const1.*dist.^2 + 2.*const2.*dist + const3)./15;
      Dd = ma2.*Dd;
      
      for i=1:length(gdist)
        ii1 = ii1+1;
        D = Dd.*gdist{i};
        D = sparse(I,J,D,n,n);
        DKff{ii1} = D + D';
      end
    else
      if isfield(gpcf, 'selectedVariables')
        x = x(:,gpcf.selectedVariables);
      end
      [n, m] =size(x);
      if ~savememory
        i1=1:m;
      else
        if i1==1
          DKff=DKff{1};
          return
        end
        i1=i1-1;
        ii1=ii1-1;
      end
      if ~isempty(gpcf.p.lengthScale)
        % loop over all the lengthScales
        if length(gpcf.lengthScale) == 1
          % In the case of isotropic PPCS3
          s2 = 1./gpcf.lengthScale.^2;
          ma2 = gpcf.magnSigma2;
          
          % Calculate the sparse distance (lower triangle) matrix
          d2 = 0;
          for i = 1:m
            d2 = d2 + s2.*(x(I,i) - x(J,i)).^2;
          end
          d = sqrt(d2);
          
          % Create the 'compact support' matrix, that is, (1-R)_+,
          % where ()_+ truncates all non-positive inputs to zero.
          cs = 1-d;
          
          % Calculate the gradient matrix
          const1 = l^3 + 9*l^2 + 23*l + 15;
          const2 = 6*l^2 + 36*l + 45;
          const3 = 15*l + 45;
          
          D = -(l+3).*cs.^(l+2).*(const1.*d.^3 + const2.*d.^2 + const3.*d + 15)/15;
          D = D + cs.^(l+3).*(3.*const1.*d.^2 + 2.*const2.*d + const3)./15;
          D = -d.*ma2.*D;
          D = sparse(I,J,D,n,n);
          
          ii1 = ii1+1;
          DKff{ii1} = D;
        else
          % In the case ARD is used
          s2 = 1./gpcf.lengthScale.^2;
          ma2 = gpcf.magnSigma2;
          
          % Calculate the sparse distance (lower triangle) matrix
          % and the distance matrix for each component
          d2 = 0;
          d_l2 = [];
          for i = 1:m
            d_l2(:,i) = s2(i).*(x(I,i) - x(J,i)).^2;
            d2 = d2 + d_l2(:,i);
          end
          d = sqrt(d2);
          d_l = d_l2;
          
          % Create the 'compact support' matrix, that is, (1-R)_+,
          % where ()_+ truncates all non-positive inputs to zero.
          cs = 1-d;
          
          const1 = l^3 + 9*l^2 + 23*l + 15;
          const2 = 6*l^2 + 36*l + 45;
          const3 = 15*l + 45;
          
          Dd = -(l+3).*cs.^(l+2).*(const1.*d.^3 + const2.*d.^2 + const3.*d + 15)/15;
          Dd = Dd + cs.^(l+3).*(3.*const1.*d.^2 + 2.*const2.*d + const3)./15;
          Dd = -ma2.*Dd;
          int = d ~= 0;
          
          for i = i1
            % Calculate the gradient matrix
            D = d_l(:,i).*Dd;
            % Divide by r in cases where r is non-zero
            D(int) = D(int)./d(int);
            D = sparse(I,J,D,n,n);
            
            ii1 = ii1+1;
            DKff{ii1} = D;
          end
        end
      end
    end
    % Evaluate the gradient of non-symmetric covariance (e.g. K_fu)
  elseif nargin == 3 || isempty(mask)
    if size(x,2) ~= size(x2,2)
      error('gpcf_ppcs -> _ghyper: The number of columns in x and x2 has to be the same. ')
    end
    
    ii1=0;
    K = gpcf.fh.cov(gpcf, x, x2);
    if ~isempty(gpcf.p.magnSigma2)
      ii1 = ii1 +1;
      DKff{ii1} = K;
    end

    l = gpcf.l;
    
    if isfield(gpcf,'metric')
      % If other than scaled euclidean metric
      [n1,m1]=size(x);
      [n2,m2]=size(x2);
      
      ma = gpcf.magnSigma2;
      
      % Compute the sparse distance matrix.
      ntriplets = nnz(K);
      I = zeros(ntriplets,1);
      J = zeros(ntriplets,1);
      R = zeros(ntriplets,1);
      dist = zeros(ntriplets,1);
      for jj = 1:length(gpcf.metric.components)
        gdist{jj} = zeros(ntriplets,1);
      end
      ntriplets = 0;
      for ii=1:n2
        d = zeros(n1,1);
        d = gpcf.metric.fh.dist(gpcf.metric, x, x2(ii,:));
        gd = gpcf.metric.fh.distg(gpcf.metric, x, x2(ii,:));
        gprior_dist = gpcf.metric.fh.lpg(gpcf.metric, x, x2(ii,:));
        
        I0t = find(d==0);
        d(d >= 1) = 0;
        [I2,J2,R2] = find(d);
        len = length(R);
        ntrip_prev = ntriplets;
        ntriplets = ntriplets + length(R2);

        ind_tr = ntrip_prev+1:ntriplets;
        I(ind_tr) = I2;
        J(ind_tr) = ii;
        dist(ind_tr) = R2;
        for jj = 1:length(gd)
          gdist{jj}(ind_tr) = gd{jj}(I2);
        end
      end

      
      ma2 = gpcf.magnSigma2;
      
      cs = 1-dist;
      
      const1 = l^3 + 9*l^2 + 23*l + 15;
      const2 = 6*l^2 + 36*l + 45;
      const3 = 15*l + 45;
      
      Dd = -(l+3).*cs.^(l+2).*(const1.*d.^3 + const2.*d.^2 + const3.*d + 15)/15;
      Dd = Dd + cs.^(l+3).*(3.*const1.*d.^2 + 2.*const2.*d + const3)./15;
      Dd = ma2.*Dd;
      
      for i=1:length(gdist)
        ii1 = ii1+1;
        D = Dd.*gdist{i};
        D = sparse(I,J,D,n1,n2);
        DKff{ii1} = D;
      end

    else
      if isfield(gpcf, 'selectedVariables')
        x = x(:,gpcf.selectedVariables);
        x2 = x2(:,gpcf.selectedVariables);
      end
      [n, m] =size(x);
      if ~savememory
        i1=1:m;
      else
        if i1==1
          DKff=DKff{1};
          return
        end
        i1=i1-1;
        ii1=ii1-1;
      end
      if ~isempty(gpcf.p.lengthScale)
        % loop over all the lengthScales
        if length(gpcf.lengthScale) == 1
          % In the case of isotropic PPCS3
          s2 = 1./gpcf.lengthScale.^2;
          ma2 = gpcf.magnSigma2;
          
          % Calculate the sparse distance (lower triangle) matrix
          dist1 = 0;
          for i=1:m
            dist1 = dist1 + s2.*(bsxfun(@minus,x(:,i),x2(:,i)')).^2;
          end
          d1 = sqrt(dist1); 
          cs1 = max(1-d1,0);
          
          const1 = l^3 + 9*l^2 + 23*l + 15;
          const2 = 6*l^2 + 36*l + 45;
          const3 = 15*l + 45;
          
          D = -(l+3).*cs1.^(l+2).*(const1.*d1.^3 + const2.*d1.^2 + const3.*d1 + 15)/15;
          D = D + cs1.^(l+3).*(3.*const1.*d1.^2 + 2.*const2.*d1 + const3)./15;
          DK_l = -d1.*ma2.*D;
          
          ii1=ii1+1;
          DKff{ii1} = DK_l;
        else
          % In the case ARD is used
          s2 = 1./gpcf.lengthScale.^2;
          ma2 = gpcf.magnSigma2;
          
          % Calculate the sparse distance (lower triangle) matrix
          % and the distance matrix for each component
          dist1 = 0; 
          d_l1 = [];
          for i = 1:m
            dist1 = dist1 + s2(i).*bsxfun(@minus,x(:,i),x2(:,i)').^2;
            d_l1{i} = s2(i).*(bsxfun(@minus,x(:,i),x2(:,i)')).^2;
          end
          d1 = sqrt(dist1); 
          cs1 = max(1-d1,0);
          
          const1 = l^3 + 9*l^2 + 23*l + 15;
          const2 = 6*l^2 + 36*l + 45;
          const3 = 15*l + 45;
          
          D = -(l+3).*cs1.^(l+2).*(const1.*d1.^3 + const2.*d1.^2 + const3.*d1 + 15)/15;
          D = D + cs1.^(l+3).*(3.*const1.*d1.^2 + 2.*const2.*d1 + const3)./15;
          
          for i = i1
            % Calculate the gradient matrix
            DK_l = -D.*ma2.*d_l1{i};
            % Divide by r in cases where r is non-zero
            DK_l(d1 ~= 0) = DK_l(d1 ~= 0)./d1(d1 ~= 0);
            ii1=ii1+1;
            DKff{ii1} = DK_l;
          end
        end
      end
    end
    % Evaluate: DKff{1}    = d mask(Kff,I) / d magnSigma2
    %           DKff{2...} = d mask(Kff,I) / d lengthScale
  elseif nargin == 4 || nargin == 5
    ii1=0;
    [n, m] =size(x);
    
    if ~isempty(gpcf.p.magnSigma2) && (~savememory || all(i1==1))
      ii1 = ii1+1;
      DKff{ii1} = gpcf.fh.trvar(gpcf, x);   % d mask(Kff,I) / d magnSigma2
    end

    if isfield(gpcf,'metric')
      dist = 0;
      distg = gpcf.metric.fh.distg(gpcf.metric, x, [], 1);
      gprior_dist = gpcf.metric.fh.lpg(gpcf.metric);
      for i=1:length(distg)
        ii1 = ii1+1;
        DKff{ii1} = 0;
      end
    else
      if ~isempty(gpcf.p.lengthScale)
        for i2=1:length(gpcf.lengthScale)
          ii1 = ii1+1;
          DKff{ii1}  = 0; % d mask(Kff,I) / d lengthScale
        end
      end
    end
  end
  if savememory
    DKff=DKff{1};
  end
end

function DKff = gpcf_ppcs3_ginput(gpcf, x, x2, i1)
%GPCF_PPCS3_GINPUT  Evaluate gradient of covariance function with 
%                   respect to x
%
%  Description
%    DKff = GPCF_PPCS3_GINPUT(GPCF, X) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_PPCS3_GINPUT(GPCF, X, X2) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2) with respect to X (cell array with matrix elements).
%    This subfunction is needed when computing gradients with 
%    respect to inducing inputs in sparse approximations.
%
%    DKff = GPCF_PPCS3_GINPUT(GPCF, X, X2, i) takes a covariance
%    function structure GPCF, a matrix X of input vectors and
%    returns DKff, the gradients of covariance matrix Kff =
%    k(X,X2), or k(X,X) if X2 is empty, with respect to ith
%    covariate in X (cell array with matrix elements). This 
%    subfunction is needed when using memory save option in
%    gp_set.
% 
%  See also
%    GPCF_PPCS3_PAK, GPCF_PPCS3_UNPAK, GPCF_PPCS3_LP, GP_G

  [n, m] =size(x);
  ii1=0;
  DKff = {};
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
  % evaluate the gradient for training covariance
  if nargin == 2 || isempty(x2)

    K = gpcf_ppcs3_trcov(gpcf, x);
    l = gpcf.l;
    [I,J] = find(K);
    
    if isfield(gpcf,'metric')
      % Compute the sparse distance matrix and its gradient.
      ntriplets = (nnz(Cdm)-n)./2;
      I = zeros(ntriplets,1);
      J = zeros(ntriplets,1);
      dist = zeros(ntriplets,1);
      for jj = 1:length(gpcf.metric.components)
        gdist{jj} = zeros(ntriplets,1);
      end
      ntriplets = 0;                
      for ii=1:n-1
        col_ind = ii + find(Cdm(ii+1:n,ii));
        d = zeros(length(col_ind),1);
        d = gpcf.metric.fh.dist(gpcf.metric, x(col_ind,:), x(ii,:));
        
        [gd, gprior_dist] = gpcf.metric.fh.ginput(gpcf.metric, x(col_ind,:), x(ii,:));

        ntrip_prev = ntriplets;
        ntriplets = ntriplets + length(d);
        
        ind_tr = ntrip_prev+1:ntriplets;
        I(ind_tr) = col_ind;
        J(ind_tr) = ii;
        dist(ind_tr) = d;
        for jj = 1:length(gd)
          gdist{jj}(ind_tr) = gd{jj};
        end
      end
      
      ma2 = gpcf.magnSigma2;
      
      cs = 1-dist;
      
      const1 = l^3 + 9*l^2 + 23*l + 15;
      const2 = 6*l^2 + 36*l + 45;
      const3 = 15*l + 45;
      
      Dd = -(l+3).*cs.^(l+2).*(const1.*d.^3 + const2.*d.^2 + const3.*d + 15)/15;
      Dd = Dd + cs.^(l+3).*(3.*const1.*d.^2 + 2.*const2.*d + const3)./15;
      Dd = ma2.*Dd;
      
      for i=1:length(gdist)
        ii1 = ii1+1;
        D = Dd.*gdist{i};
        D = sparse(I,J,D,n,n);
        DKff{ii1} = D + D';
      end
      
    else
      if length(gpcf.lengthScale) == 1
        % In the case of an isotropic
        s2 = repmat(1./gpcf.lengthScale.^2, 1, m);
      else
        s2 = 1./gpcf.lengthScale.^2;
      end
      ma2 = gpcf.magnSigma2;
      
      % Calculate the sparse distance (lower triangle) matrix
      % and the distance matrix for each component
      d2 = 0;
      for i = 1:m
        d2 = d2 + s2(i).*(x(I,i) - x(J,i)).^2;
      end
      d = sqrt(d2);
      
      % Create the 'compact support' matrix, that is, (1-R)_+,
      % where ()_+ truncates all non-positive inputs to zero.
      cs = 1-d;
      
      const1 = l^3 + 9*l^2 + 23*l + 15;
      const2 = 6*l^2 + 36*l + 45;
      const3 = 15*l + 45;
      
      Dd = -(l+3).*cs.^(l+2).*(const1.*d.^3 + const2.*d.^2 + const3.*d + 15)/15;
      Dd = Dd + cs.^(l+3).*(3.*const1.*d.^2 + 2.*const2.*d + const3)./15;
      Dd = ma2.*Dd;
      
      Dd = sparse(I,J,Dd,n,n);
      d = sparse(I,J,d,n,n);
      
      row = ones(n,1);
      cols = 1:n;
      for i = i1
        for j = 1:n
          % Calculate the gradient matrix
          ind = find(d(:,j));
          apu = full(Dd(:,j)).*s2(i).*(x(j,i)-x(:,i));
          apu(ind) = apu(ind)./d(ind,j);
          D = sparse(row*j, cols, apu, n, n);
          D = D+D';
          
          ii1 = ii1+1;
          DKff{ii1} = D;
        end
      end
    end

    % Evaluate the gradient of non-symmetric covariance (e.g. K_fu)
  elseif nargin == 3 || nargin == 4
    if size(x,2) ~= size(x2,2)
      error('gpcf_ppcs -> _ghyper: The number of columns in x and x2 has to be the same. ')
    end
    
    ii1=0;
    K = gpcf.fh.cov(gpcf, x, x2);
    n2 = size(x2,1);
    l = gpcf.l;
    
    if isfield(gpcf,'metric')
      % If other than scaled euclidean metric
      [n1,m1]=size(x);
      [n2,m2]=size(x2);
      
      ma = gpcf.magnSigma2;
      
      % Compute the sparse distance matrix.
      ntriplets = nnz(K);
      I = zeros(ntriplets,1);
      J = zeros(ntriplets,1);
      R = zeros(ntriplets,1);
      dist = zeros(ntriplets,1);
      for jj = 1:length(gpcf.metric.components)
        gdist{jj} = zeros(ntriplets,1);
      end
      ntriplets = 0;
      for ii=1:n2
        d = zeros(n1,1);
        d = gpcf.metric.fh.dist(gpcf.metric, x, x2(ii,:));
        [gd, gprior_dist] = gpcf.metric.fh.ginput(gpcf.metric, x, x2(ii,:));
        
        I0t = find(d==0);
        d(d >= 1) = 0;
        [I2,J2,R2] = find(d);
        len = length(R);
        ntrip_prev = ntriplets;
        ntriplets = ntriplets + length(R2);

        ind_tr = ntrip_prev+1:ntriplets;
        I(ind_tr) = I2;
        J(ind_tr) = ii;
        dist(ind_tr) = R2;
        for jj = 1:length(gd)
          gdist{jj}(ind_tr) = gd{jj}(I2);
        end
      end

      
      ma2 = gpcf.magnSigma2;
      
      cs = 1-dist;
      
      const1 = l^3 + 9*l^2 + 23*l + 15;
      const2 = 6*l^2 + 36*l + 45;
      const3 = 15*l + 45;
      
      Dd = -(l+3).*cs.^(l+2).*(const1.*d.^3 + const2.*d.^2 + const3.*d + 15)/15;
      Dd = Dd + cs.^(l+3).*(3.*const1.*d.^2 + 2.*const2.*d + const3)./15;
      Dd = ma2.*Dd;
      
      for i=1:length(gdist)
        ii1 = ii1+1;
        D = Dd.*gdist{i};
        D = sparse(I,J,D,n1,n2);
        DKff{ii1} = D;
      end

    else
      if length(gpcf.lengthScale) == 1
        % In the case of an isotropic PPCS3
        s2 = repmat(1./gpcf.lengthScale.^2, 1, m);
      else
        s2 = 1./gpcf.lengthScale.^2;
      end
      ma2 = gpcf.magnSigma2;
      
      % Calculate the sparse distance (lower triangle) matrix
      % and the distance matrix for each component
      d2 = 0; 
      for i = 1:m
        d2 = d2 + s2(i).*bsxfun(@minus,x(:,i),x2(:,i)').^2;
      end
      d = sqrt(d2);
      cs = max(1-d,0);

      const1 = l^3 + 9*l^2 + 23*l + 15;
      const2 = 6*l^2 + 36*l + 45;
      const3 = 15*l + 45;
      
      Dd = -(l+3).*cs.^(l+2).*(const1.*d.^3 + const2.*d.^2 + const3.*d + 15)/15;
      Dd = Dd + cs.^(l+3).*(3.*const1.*d.^2 + 2.*const2.*d + const3)./15;
      Dd = ma2.*Dd;

      row = ones(n2,1);
      cols = 1:n2;
      for i = i1
        for j = 1:n
          % Calculate the gradient matrix
          ind = find(d(j,:));
          apu = Dd(j,:).*s2(i).*(x(j,i)-x2(:,i))';
          apu(ind) = apu(ind)./d(j,ind);
          D = sparse(row*j, cols, apu, n, n2);
          
          ii1 = ii1+1;
          DKff{ii1} = D;
        end
      end
    end
  end
end


function C = gpcf_ppcs3_cov(gpcf, x1, x2, varargin)
%GP_PPCS3_COV  Evaluate covariance matrix between two input vectors
%
%  Description         
%    C = GP_PPCS3_COV(GP, TX, X) takes in covariance function of
%    a Gaussian process GP and two matrixes TX and X that contain
%    input vectors to GP. Returns covariance matrix C. Every
%    element ij of C contains covariance between inputs i in TX
%    and j in X. This is a mandatory subfunction used for example in
%    prediction and energy computations.
%
%  See also
%    GPCF_PPCS3_TRCOV, GPCF_PPCS3_TRVAR, GP_COV, GP_TRCOV

  if isfield(gpcf,'metric')
    % If other than scaled euclidean metric
    [n1,m1]=size(x1);
    [n2,m2]=size(x2);
  else
    % If scaled euclidean metric
    if isfield(gpcf, 'selectedVariables')
      x1 = x1(:,gpcf.selectedVariables);
      x2 = x2(:,gpcf.selectedVariables);
    end
    [n1,m1]=size(x1);
    [n2,m2]=size(x2);
    
    s = 1./(gpcf.lengthScale);
    s2 = s.^2;
    if size(s)==1
      s2 = repmat(s2,1,m1);
    end
  end
  ma2 = gpcf.magnSigma2;
  l = gpcf.l;
  
  % Compute the sparse distance matrix.
  ntriplets = max(1,floor(0.03*n1*n2));
  I = zeros(ntriplets,1);
  J = zeros(ntriplets,1);
  R = zeros(ntriplets,1);
  ntriplets = 0;
  I0=zeros(ntriplets,1);
  J0=zeros(ntriplets,1);
  nn0=0;
  for ii1=1:n2
    d = zeros(n1,1);
    if isfield(gpcf, 'metric')
      d = gpcf.metric.fh.dist(gpcf.metric, x1, x2(ii1,:));
    else
      for j=1:m1
        d = d + s2(j).*(x1(:,j)-x2(ii1,j)).^2;
      end
    end
    %d = sqrt(d);
    I0t = find(d==0);
    d(d >= 1) = 0;
    [I2,J2,R2] = find(d);
    R2=sqrt(R2);
    %len = length(R);
    ntrip_prev = ntriplets;
    ntriplets = ntriplets + length(R2);
    
    I(ntrip_prev+1:ntriplets) = I2;
    J(ntrip_prev+1:ntriplets) = ii1;
    R(ntrip_prev+1:ntriplets) = R2;
    I0(nn0+1:nn0+length(I0t)) = I0t;
    J0(nn0+1:nn0+length(I0t)) = ii1;
    nn0 = nn0+length(I0t);
  end
  r = sparse(I(1:ntriplets),J(1:ntriplets),R(1:ntriplets));
  [I,J,r] = find(r);
  cs = full(sparse(max(0, 1-r)));
  
  const1 = l^3 + 9*l^2 + 23*l + 15;
  const2 = 6*l^2 + 36*l + 45;
  const3 = 15*l + 45;
  C = ma2.*cs.^(l+3).*(const1.*r.^3 + const2.*r.^2 + const3.*r + 15)/15;
  
  C = sparse(I,J,C,n1,n2) + sparse(I0,J0,ma2,n1,n2);
end

function C = gpcf_ppcs3_trcov(gpcf, x)
%GP_PPCS3_TRCOV  Evaluate training covariance matrix of inputs
%
%  Description
%    C = GP_PPCS3_TRCOV(GP, TX) takes in covariance function of a
%    Gaussian process GP and matrix TX that contains training
%    input vectors. Returns covariance matrix C. Every element ij
%    of C contains covariance between inputs i and j in TX. This is 
%    a mandatory subfunction used for example in prediction and 
%    energy computations.
%
%  See also
%    GPCF_PPCS3_COV, GPCF_PPCS3_TRVAR, GP_COV, GP_TRCOV

  if isfield(gpcf,'metric') 
    % If other than scaled euclidean metric
    [n, m] =size(x);            
  else
    % If a scaled euclidean metric try first mex-implementation 
    % and if there is not such... 
    C = trcov(gpcf,x);
    % ... evaluate the covariance here.
    if isnan(C)
      if isfield(gpcf, 'selectedVariables')
        x = x(:,gpcf.selectedVariables);
      end
      [n, m] =size(x);
      
      s = 1./(gpcf.lengthScale);
      s2 = s.^2;
      if size(s)==1
        s2 = repmat(s2,1,m);
      end
    else
      return
    end
  end
  ma2 = gpcf.magnSigma2;
  l = gpcf.l;
  
  % Compute the sparse distance matrix.
  ntriplets = max(1,floor(0.03*n*n));
  I = zeros(ntriplets,1);
  J = zeros(ntriplets,1);
  R = zeros(ntriplets,1);
  ntriplets = 0;
  ntripletsz = max(1,floor(0.03.^2*n*n));
  Iz = zeros(ntripletsz,1);
  Jz = zeros(ntripletsz,1);
  ntripletsz = 0;
  for ii1=1:n-1
    d = zeros(n-ii1,1);
    col_ind = ii1+1:n;
    if isfield(gpcf,'metric')
      d = gpcf.metric.fh.dist(gpcf.metric, x(col_ind,:), x(ii1,:));
    else
      for ii2=1:m
        d = d+s2(ii2).*(x(col_ind,ii2)-x(ii1,ii2)).^2;
      end
    end
    %d = sqrt(d);
    % store zero distance index
    [I2z,J2z] = find(d==0);
    % create triplets for distances 0<d<1
    d(d >= 1) = 0;
    [I2,J2,R2] = find(d);
    len = length(R);
    ntrip_prev = ntriplets;
    ntriplets = ntriplets + length(R2);
    if (ntriplets > len)
      I(2*len) = 0;
      J(2*len) = 0;
      R(2*len) = 0;
    end
    ind_tr = ntrip_prev+1:ntriplets;
    I(ind_tr) = ii1+I2;
    J(ind_tr) = ii1;
    R(ind_tr) = sqrt(R2);
    % create triplets for distances d==0 (i~=j)
    lenz = length(Iz);
    ntrip_prevz = ntripletsz;
    ntripletsz = ntripletsz + length(I2z);
    if (ntripletsz > lenz)
      Iz(2*lenz) = 0;
      Jz(2*lenz) = 0;
    end
    ind_trz = ntrip_prevz+1:ntripletsz;
    Iz(ind_trz) = ii1+I2z;
    Jz(ind_trz) = ii1;
  end
  % create a lower triangular sparse distance matrix from the triplets for distances 0<d<1
  R = sparse(I(1:ntriplets),J(1:ntriplets),R(1:ntriplets),n,n);
  % create a lower triangular sparse covariance matrix from the
  % triplets for distances d==0 (i~=j)
  Rz = sparse(Iz(1:ntripletsz),Jz(1:ntripletsz),repmat(ma2,1,ntripletsz),n,n);
  
  % Find the non-zero elements of R.
  [I,J,rn] = find(R);
  % Compute covariances for distances 0<d<1
  const1 = l^3 + 9*l^2 + 23*l + 15;
  const2 = 6*l^2 + 36*l + 45;
  const3 = 15*l + 45;
  cs = max(0,1-rn);
  C = ma2.*cs.^(l+3).*(const1.*rn.^3 + const2.*rn.^2 + const3.*rn + 15)/15;
  % create a lower triangular sparse covariance matrix from the triplets for distances 0<d<1
  C = sparse(I,J,C,n,n);
  % add the lower triangular covariance matrix for distances d==0 (i~=j)
  C = C + Rz;
  % form a square covariance matrix and add the covariance matrix for i==j (d==0)
  C = C + C' + sparse(1:n,1:n,ma2,n,n);

end

function C = gpcf_ppcs3_trvar(gpcf, x)
%GP_PPCS3_TRVAR  Evaluate training variance vector
%
%  Description
%    C = GP_PPCS3_TRVAR(GPCF, TX) takes in covariance function of
%    a Gaussian process GPCF and matrix TX that contains training
%    inputs. Returns variance vector C. Every element i of C
%    contains variance of input i in TX. This is a mandatory 
%    subfunction used for example in prediction and energy computations.
%
%  See also
%    GPCF_PPCS3_COV, GP_COV, GP_TRCOV

  [n, m] =size(x);

  C = ones(n,1).*gpcf.magnSigma2;
  C(C<eps)=0;
end

function reccf = gpcf_ppcs3_recappend(reccf, ri, gpcf)
%RECAPPEND  Record append
%
%  Description
%    RECCF = GPCF_PPCS3_RECAPPEND(RECCF, RI, GPCF) takes a
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
    reccf.type = 'gpcf_ppcs3';
    reccf.nin = ri.nin;
    reccf.l = floor(reccf.nin/2)+4;

    % cf is compactly supported
    reccf.cs = 1;
    
    % Initialize parameters
    reccf.lengthScale= [];
    reccf.magnSigma2 = [];

    % Set the function handles
    reccf.fh.pak = @gpcf_ppcs3_pak;
    reccf.fh.unpak = @gpcf_ppcs3_unpak;
    reccf.fh.e = @gpcf_ppcs3_lp;
    reccf.fh.lpg = @gpcf_ppcs3_lpg;
    reccf.fh.cfg = @gpcf_ppcs3_cfg;
    reccf.fh.cov = @gpcf_ppcs3_cov;
    reccf.fh.trcov  = @gpcf_ppcs3_trcov;
    reccf.fh.trvar  = @gpcf_ppcs3_trvar;
    reccf.fh.recappend = @gpcf_ppcs3_recappend;  
    reccf.p=[];
    reccf.p.lengthScale=[];
    reccf.p.magnSigma2=[];
    if isfield(ri.p,'lengthScale') && ~isempty(ri.p.lengthScale)
      reccf.p.lengthScale = ri.p.lengthScale;
    end
    if ~isempty(ri.p.magnSigma2)
      reccf.p.magnSigma2 = ri.p.magnSigma2;
    end
    if isfield(ri, 'selectedVariables')
        reccf.selectedVariables = ri.selectedVariables;
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
  
  end
end

