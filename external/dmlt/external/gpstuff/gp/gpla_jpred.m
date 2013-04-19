function [Eft, Covft, ljpyt] = gpla_jpred(gp, x, y, varargin)
%GPLA_PRED  Predictions with Gaussian Process Laplace approximation
%
%  Description
%    [EFT, COVFT] = GPLA_JPRED(GP, X, Y, XT, OPTIONS)
%    takes a GP structure together with matrix X of training
%    inputs and vector Y of training targets, and evaluates the
%    predictive distribution at test inputs XT. Returns a posterior
%    mean EFT and covariance COVFT of latent variables.
%
%        Eft =  E[f | xt,x,y,th]  = K_fy*(Kyy+s^2I)^(-1)*y
%      Covft = Cov[f | xt,x,y,th] = K_fy - K_fy*(Kyy+s^2I)^(-1)*K_yf. 
%
%    Each row of X corresponds to one input vector and each row of
%    Y corresponds to one output vector.
%
%    [EFT, COVFT, LJPYT] = GPLA_JPRED(GP, X, Y, XT, 'yt', YT, ...) 
%    returns also logarithm of the predictive joint density JPYT of
%    the observations YT at test input locations XT. This can be
%    used for example in the cross-validation. Here Y has to be
%    vector.
%
%    [EFT, COVFT, LJPYT, EYT, VARYT] = GPLA_JPRED(GP, X, Y, XT, 'yt', YT, ...)
%    returns also the posterior predictive mean and covariance.
%
%    [EF, COVF, LJPY] = GPEP_JPRED(GP, X, Y, OPTIONS)
%    evaluates the predictive distribution at training inputs X
%    and logarithm of the predictive density PY of the training
%    observations Y.
%
%    OPTIONS is optional parameter-value pair
%      predcf - an index vector telling which covariance functions are 
%               used for prediction. Default is all (1:gpcfn). 
%               See additional information below.
%      tstind - a vector/cell array defining, which rows of X belong 
%               to which training block in *IC type sparse models. 
%               Default is []. In case of PIC, a cell array
%               containing index vectors specifying the blocking
%               structure for test data. IN FIC and CS+FIC a
%               vector of length n that points out the test inputs
%               that are also in the training set (if none, set
%               TSTIND = [])
%      yt     - optional observed yt in test points (see below)
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case. 
%      zt     - optional observed quantity in triplet (xt_i,yt_i,zt_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, the expected 
%               value for the ith case. 
%
%    NOTE! In case of FIC and PIC sparse approximation the
%    prediction for only some PREDCF covariance functions is just
%    an approximation since the covariance functions are coupled in
%    the approximation and are not strictly speaking additive
%    anymore.
%
%    For example, if you use covariance such as K = K1 + K2 your
%    predictions Eft1 = ep_pred(GP, X, Y, X, 'predcf', 1) and 
%    Eft2 = ep_pred(gp, x, y, x, 'predcf', 2) should sum up to 
%    Eft = ep_pred(gp, x, y, x). That is Eft = Eft1 + Eft2. With 
%    FULL model this is true but with FIC and PIC this is true only 
%    approximately. That is Eft \approx Eft1 + Eft2.
%
%    With CS+FIC the predictions are exact if the PREDCF covariance
%    functions are all in the FIC part or if they are CS
%    covariances.
%
%    NOTE! When making predictions with a subset of covariance
%    functions with FIC approximation the predictive variance can
%    in some cases be ill-behaved i.e. negative or unrealistically
%    small. This may happen because of the approximative nature of
%    the prediction.
%
%  See also
%    GPLA_E, GPLA_G, GP_PRED, DEMO_SPATIAL, DEMO_CLASSIFIC

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2011-2012 Ville Tolvanen
% Copyright (c) 2012 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPLA_JPRED';
  ip.addRequired('gp', @isstruct);
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addOptional('xt', [], @(x) isempty(x) || (isreal(x) && all(isfinite(x(:)))))
  ip.addParamValue('yt', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('zt', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('predcf', [], @(x) isempty(x) || ...
                   isvector(x) && isreal(x) && all(isfinite(x)&x>0))
  ip.addParamValue('tstind', [], @(x) isempty(x) || iscell(x) ||...
                   (isvector(x) && isreal(x) && all(isfinite(x)&x>0)))
  if numel(varargin)==0 || isnumeric(varargin{1})
    % inputParser should handle this, but it doesn't
    ip.parse(gp, x, y, varargin{:});
  else
    ip.parse(gp, x, y, [], varargin{:});
  end
  xt=ip.Results.xt;
  yt=ip.Results.yt;
  z=ip.Results.z;
  zt=ip.Results.zt;
  predcf=ip.Results.predcf;
  tstind=ip.Results.tstind;
  if isempty(xt)
    xt=x;
    if isempty(tstind)
      if iscell(gp)
        gptype=gp{1}.type;
      else
        gptype=gp.type;
      end
      switch gptype
        case {'FULL' 'VAR' 'DTC' 'SOR'}
          tstind = [];
        case {'FIC' 'CS+FIC'}
          tstind = 1:size(x,1);
        case 'PIC'
          if iscell(gp)
            tstind = gp{1}.tr_index;
          else
            tstind = gp.tr_index;
          end
      end
    end
    if isempty(yt)
      yt=y;
    end
    if isempty(zt)
      zt=z;
    end
  end

  [tn, tnin] = size(x);
  
  switch gp.type
    case 'FULL'
      % ============================================================
      % FULL
      % ============================================================
      [e, edata, eprior, f, L, a, W, p] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);

      ntest=size(xt,1);
      K_nf = gp_cov(gp,xt,x,predcf);
      if isfield(gp,'meanf')
        [H,b_m,B_m Hs]=mean_prep(gp,x,xt);
        K_nf=K_nf + Hs'*B_m*H;
        K = gp_trcov(gp, x);
        K = K+H'*B_m*H;
      end
      
      % Evaluate the mean
      if issparse(K_nf) && issparse(L)        
        deriv = gp.lik.fh.llg(gp.lik, y(p), f, 'latent', z(p));
        Eft = K_nf(:,p)*deriv;
      else
        deriv = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
        Eft = K_nf*deriv;
        if isfield(gp,'meanf')
          Eft=Eft + K_nf*(K\Hs'*b_m);
        end
      end

      if nargout > 1
        % Evaluate the variance
        %           kstarstar = gp_trvar(gp,xt,predcf);
        kstarstar = gp_trcov(gp,xt,predcf);               
        if isfield(gp,'meanf')
          kstarstar= kstarstar + diag(Hs'*B_m*Hs);
        end
        if W >= 0
          % This is the usual case where likelihood is log concave
          % for example, Poisson and probit
          if issparse(K_nf) && issparse(L)
            % If compact support covariance functions are used 
            % the covariance matrix will be sparse
            K = gp_trcov(gp, x);                            
            sqrtW = sparse(1:tn, 1:tn, sqrt(W), tn, tn);
            sqrtWKfn = sqrtW*K_nf';
            V = ldlsolve(L,sqrtWKfn);
            Covft = K - sqrtWKfn'*V;
          else
            W = diag(W);
            V = L\(sqrt(W)*K_nf');
            Covft = kstarstar - (V'*V);
          end
        else                  
          % We may end up here if the likelihood is not log concace
          % For example Student-t likelihood
          V = L*diag(W);
          R = diag(W) - V'*V;
          Covft = kstarstar - K_nf*(R*K_nf');
        end
      end
      
    case 'FIC'        
      % ============================================================
      % FIC
      % ============================================================    
      % Predictions with FIC sparse approximation for GP
      % Here tstind = 1 if the prediction is made for the training set 
      if nargin > 6
        if ~isempty(tstind) && length(tstind) ~= size(x,1)
          error('tstind (if provided) has to be of same lenght as x.')
        end
      else
        tstind = [];
      end

      u = gp.X_u;
      K_fu = gp_cov(gp, x, u, predcf);         % f x u
      K_uu = gp_trcov(gp, u, predcf);          % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;                  % ensure the symmetry of K_uu
      Luu = chol(K_uu)';

      m = size(u,1);

      [e, edata, eprior, f, L, a, La2] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);

      deriv = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
      ntest=size(xt,1);

      K_nu=gp_cov(gp,xt,u,predcf);
      Eft = K_nu*(Luu'\(Luu\(K_fu'*deriv)));

      % if the prediction is made for training set, evaluate Lav also for prediction points
      if ~isempty(tstind)
        [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf);
        B=Luu\(K_fu');
        Qv_ff=sum(B.^2)';
        %Lav = zeros(size(La));
        %Lav(tstind) = Kv_ff-Qv_ff;
        Lav = Kv_ff-Qv_ff;            
        Eft(tstind) = Eft(tstind) + Lav.*deriv;
      end

      
      % Evaluate the variance
      if nargout > 1
        % re-evaluate matrices with training components
        Kfu_tr = gp_cov(gp, x, u);
        Kuu_tr = gp_trcov(gp, u);
        Kuu_tr = (K_uu+K_uu')./2;
        
        W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
        kstarstar = gp_trvar(gp,xt,predcf);
        La = W.*La2;
        Lahat = 1 + La;
        B = (repmat(sqrt(W),1,m).*Kfu_tr);

        % Components for (I + W^(1/2)*(Qff + La2)*W^(1/2))^(-1) = Lahat^(-1) - L2*L2'
        B2 = repmat(Lahat,1,m).\B;
        A2 = Kuu_tr + B'*B2; A2=(A2+A2)/2;
        L2 = B2/chol(A2);
        
        
        % NOTE!                                             % MUUTETTU
        % This is done with full matrices at the moment. 
        % Needs to be rewritten.
        if isempty(tstind)
          error('tstind not provided.');
        end
        [Knn_v, Cnn_v] = gp_trvar(gp,xt,predcf);
        Luu = chol(K_uu)';
        B=Luu\(K_fu');
        B2 = Luu\(K_nu');
        Lav_n = Knn_v - sum(B2.^2)';
        BL = B*L;
        
        K_nf = B2'*B + diag(Lav);
        K = B2'*B2 + diag(Lav_n);
        C = B'*B + diag(Lav);
        W = diag(W);
        B = eye(size(C)) + sqrt(W)*C*sqrt(W);
        L = chol(B)';
        
        V = L\(sqrt(W)*K_nf');
        Covft = K - V'*V;

        %             % Set params for K_nf
        %             BB=Luu\(B');
        %             BB2=Luu\(K_nu');
        %             Covft = kstarstar - sum(BB2'.*(BB*(repmat(Lahat,1,m).\BB')*BB2)',2)  + sum((K_nu*(K_uu\(B'*L2))).^2, 2);
        %             
        %             % if the prediction is made for training set, evaluate Lav also for prediction points
        %             if ~isempty(tstind)
        %                 LavsW = Lav.*sqrt(W);
        %                     Covft(tstind) = Covft(tstind) - (LavsW./sqrt(Lahat)).^2 + sum((repmat(LavsW,1,m).*L2).^2, 2) ...
        %                            - 2.*sum((repmat(LavsW,1,m).*(repmat(Lahat,1,m).\B)).*(K_uu\K_nu(tstind,:)')',2)...
        %                            + 2.*sum((repmat(LavsW,1,m).*L2).*(L2'*B*(K_uu\K_nu(tstind,:)'))' ,2);
        %             end
      end
      
    case {'PIC' 'PIC_BLOCK'}        
      % ============================================================
      % PIC
      % ============================================================
      % Predictions with PIC sparse approximation for GP
      u = gp.X_u;
      K_fu = gp_cov(gp, x, u, predcf);         % f x u
      K_uu = gp_trcov(gp, u, predcf);          % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;                  % ensure the symmetry of K_uu
      K_nu=gp_cov(gp,xt,u,predcf);

      ind = gp.tr_index;
      ntest = size(xt,1);
      m = size(u,1);
      Luu = chol(K_uu)';
      B=Luu\(K_fu');
      B2 = Luu\(K_nu');        

      [e, edata, eprior, f, L, a, La2] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);

      deriv = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);

      iKuuKuf = K_uu\K_fu';
      w_bu=zeros(length(xt),length(u));
      w_n=zeros(length(xt),1);
      for i=1:length(ind)
        w_bu(tstind{i},:) = repmat((iKuuKuf(:,ind{i})*deriv(ind{i},:))', length(tstind{i}),1);
        K_nf = gp_cov(gp, xt(tstind{i},:), x(ind{i},:), predcf);              % n x u
        w_n(tstind{i},:) = K_nf*deriv(ind{i},:);
      end

      Eft = K_nu*(iKuuKuf*deriv) - sum(K_nu.*w_bu,2) + w_n;

      % Evaluate the variance
      if nargout > 1
        W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
        kstarstar = gp_trvar(gp,xt,predcf);
        sqrtW = sqrt(W);
        % Components for (I + W^(1/2)*(Qff + La2)*W^(1/2))^(-1) = Lahat^(-1) - L2*L2'
        for i=1:length(ind)
          La{i} = diag(sqrtW(ind{i}))*La2{i}*diag(sqrtW(ind{i}));
          Lahat{i} = eye(size(La{i})) + La{i};
        end
        %             B = (repmat(sqrt(W),1,m).*K_fu);
        sKfu = (repmat(sqrt(W),1,m).*K_fu);
        for i=1:length(ind)
          iLasKfu(ind{i},:) = Lahat{i}\sKfu(ind{i},:);
        end
        A2 = K_uu + sKfu'*iLasKfu; A2=(A2+A2)/2;
        L2 = iLasKfu/chol(A2);

        % NOTE!                                             
        % This is done with full matrices at the moment. 
        % Needs to be rewritten.
        Knn = B2'*B2;
        Knf = B2'*B;
        C = -L2*L2';
        for i=1:length(ind)
          La = gp_trcov(gp, xt(tstind{i},:), predcf) - B2(:,tstind{i})'*B2(:,tstind{i});
          Knn(ind{i},ind{i}) =  Knn(ind{i},ind{i}) + La;
          Laa = gp_cov(gp, xt(tstind{i},:), x(ind{i},:),predcf) - B2(:,tstind{i})'*B(:,ind{i});
          Knf(tstind{i},ind{i}) =  Knf(tstind{i},ind{i}) + Laa;
          C(ind{i},ind{i}) =  C(ind{i},ind{i}) + inv(Lahat{i});
        end
        C = diag(sqrtW)*C*diag(sqrtW);
        
        Covft = Knn - Knf * C * Knf';
      end
      
    case 'CS+FIC'        
      % ============================================================
      % CS+FIC
      % ============================================================
      % Predictions with CS+FIC sparse approximation for GP
      % Here tstind = 1 if the prediction is made for the training set 
      if nargin > 6
        if ~isempty(tstind) && length(tstind) ~= size(x,1)
          error('tstind (if provided) has to be of same lenght as x.')
        end
        %         else
        %              tstind = [];
      end

      n = size(x,1);
      n2 = size(xt,1);
      u = gp.X_u;
      m = length(u);

      [e, edata, eprior, f, L, a, La2] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
      
      % Indexes to all non-compact support and compact support covariances.
      cf1 = [];
      cf2 = [];
      % Indexes to non-CS and CS covariances, which are used for predictions
      predcf1 = [];
      predcf2 = [];    
      
      ncf = length(gp.cf);
      % Loop through all covariance functions
      for i = 1:ncf        
        % Non-CS covariances
        if ~isfield(gp.cf{i},'cs') 
          cf1 = [cf1 i];
          % If used for prediction
          if ~isempty(find(predcf==i))
            predcf1 = [predcf1 i]; 
          end
          % CS-covariances
        else
          cf2 = [cf2 i];           
          % If used for prediction
          if ~isempty(find(predcf==i))
            predcf2 = [predcf2 i]; 
          end
        end
      end
      if isempty(predcf1) && isempty(predcf2)
        predcf1 = cf1;
        predcf2 = cf2;
      end
      
      % Determine the types of the covariance functions used
      % in making the prediction.
      if ~isempty(predcf1) && isempty(predcf2)       % Only non-CS covariances
        ptype = 1;
        predcf2 = cf2;
      elseif isempty(predcf1) && ~isempty(predcf2)   % Only CS covariances
        ptype = 2;
        predcf1 = cf1;
      else                                           % Both non-CS and CS covariances
        ptype = 3;
      end
      
      K_fu = gp_cov(gp,x,u,predcf1);          % f x u
      K_uu = gp_trcov(gp,u,predcf1);          % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;                 % ensure the symmetry of K_uu
      K_nu=gp_cov(gp,xt,u,predcf1);

      Kcs_nf = gp_cov(gp, xt, x, predcf2);
      Kcs_nn = gp_trcov(gp, xt, predcf2);

      deriv = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
      ntest=size(xt,1);

      % Calculate the predictive mean according to the type of
      % covariance functions used for making the prediction
      if ptype == 1
        Eft = K_nu*(K_uu\(K_fu'*deriv));
      elseif ptype == 2
        Eft = Kcs_nf*deriv;
      else 
        Eft = K_nu*(K_uu\(K_fu'*deriv)) + Kcs_nf*deriv;        
      end

      % evaluate also Lav if the prediction is made for training set
      if ~isempty(tstind)
        [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf1);
        Luu = chol(K_uu)';
        B=Luu\(K_fu');
        Qv_ff=sum(B.^2)';
        %Lav = zeros(size(Eft));
        %Lav(tstind) = Kv_ff-Qv_ff;
        Lav = Kv_ff-Qv_ff;
      end

      % Add also Lav if the prediction is made for training set
      % and non-CS covariance function is used for prediction
      if ~isempty(tstind) && (ptype == 1 || ptype == 3)
        Eft(tstind) = Eft(tstind) + Lav.*deriv;
      end
      
      % Evaluate the variance
      if nargout > 1
        W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
        sqrtW = sparse(1:tn,1:tn,sqrt(W),tn,tn);
        kstarstar = gp_trvar(gp,xt,predcf);
        Luu = chol(K_uu)';
        Lahat = sparse(1:tn,1:tn,1,tn,tn) + sqrtW*La2*sqrtW;
        B = sqrtW*K_fu;

        % Components for (I + W^(1/2)*(Qff + La2)*W^(1/2))^(-1) = Lahat^(-1) - L2*L2'
        B2 = Lahat\B;
        A2 = K_uu + B'*B2; A2=(A2+A2)/2;
        L2 = B2/chol(A2);

        % Set params for K_nf
        BB=Luu\(B)';    % sqrtW*K_fu
        BB2=Luu\(K_nu');
        iLaKfu = La2\K_fu;
        
        [Knn_v, Cnn_v] = gp_trvar(gp,xt,predcf1);
        % Calculate the predictive variance according to the type
        % of covariance functions used for making the prediction
        
        % NOTE! We form full matrices below. Needs to be rewritten so that at most nxm matrices are formed
        
        % Check results with full matrices    
        Lav_pr = Kv_ff-Qv_ff;
        B=Luu\(K_fu');
        B2=Luu\(K_nu');
        
        K = B2'*B2 + Kcs_nn + diag(Knn_v - sum(B2.^2)');
        C = B'*B + Kcs_nn + diag(Lav);
        K_nf = B2'*B + Kcs_nf + diag(Lav_pr); 
        
        W = diag(W);
        B = eye(size(C)) + sqrt(W)*C*sqrt(W);
        L = chol(B)';
        
        V = L\(sqrt(W)*K_nf');
        Covft = K - V'*V;

      end
  end
  
  if nargout > 2
    [sampft] = gp_rnd(gp,x,y, xt, 'z', z, 'zt', zt, 'nsamp', 500);
    lpyt = zeros(500,1);
    for i=1:size(sampft,2)
      lpyt(i) = gp.lik.fh.ll(gp.lik, y, sampft(:,i), z);
    end
    ljpyt = sumlogs(lpyt);
  end
  
  if nargout > 3
    error('Too many output arguments for GPLA_JPRED.')
  end
end
