function [Eft, Covft, ljpyt, Eyt, Covyt] = gp_jpred(gp, x, y, varargin)
%GP_JPRED  Joint predictions with Gaussian process 
%
%  Description
%    [EFT, COVFT] = GP_JPRED(GP, X, Y, XT, OPTIONS)
%    takes a GP structure together with matrix X of training
%    inputs and vector Y of training targets, and evaluates the
%    joint predictive distribution at test inputs XT. Returns a posterior
%    mean EFT and covariance COVFT of latent variables.
%
%        Eft =  E[f | xt,x,y,th]  = K_fy*(Kyy+s^2I)^(-1)*y
%      Covft = Var[f | xt,x,y,th] = K_fy - K_fy*(Kyy+s^2I)^(-1)*K_yf. 
%
%    Each row of X corresponds to one input vector and each row of
%    Y corresponds to one output vector.
%
%    [EFT, COVFT, LJPYT] = GP_JPRED(GP, X, Y, XT, 'yt', YT, ...)
%    returns also logarithm of the predictive joint density LJPYT of
%    the observations YT at test input locations XT. This can be
%    used for example in the cross-validation. Here Y has to be
%    vector.
%
%    [EFT, COVFT, LJPYT, EYT, COVYT] = GP_JPRED(GP, X, Y, XT, 'yt', YT, ...)
%    returns also the posterior predictive mean and covariance.
% 
%    [EF, COVF, LJPY, EY, VARY] = GP_JPRED(GP, X, Y, OPTIONS)
%    evaluates the predictive distribution at training inputs X
%    and logarithm of the predictive density PY of the training
%    observations Y.
%
%    OPTIONS is optional parameter-value pair
%      predcf - an index vector telling which covariance functions are 
%                 used for prediction. Default is all (1:gpcfn). 
%                 See additional information below.
%      tstind - a vector/cell array defining, which rows of X belong 
%               to which training block in *IC type sparse models. 
%               Default is []. In case of PIC, a cell array
%               containing index vectors specifying the blocking
%               structure for test data. IN FIC and CS+FIC a
%               vector of length n that points out the test inputs
%               that are also in the training set (if none, set
%               TSTIND = [])
%      yt     - optional observed yt in test points (see below)
%
%    NOTE! In case of FIC and PIC sparse approximation the
%    prediction for only some PREDCF covariance functions is just
%    an approximation since the covariance functions are coupled in
%    the approximation and are not strictly speaking additive
%    anymore.
%
%    For example, if you use covariance such as K = K1 + K2 your
%    predictions Ef1 = gp_pred(GP, X, Y, X, 'predcf', 1) and Ef2 =
%    gp_pred(gp, x, y, x, 'predcf', 2) should sum up to Ef =
%    gp_pred(gp, x, y, x). That is Ef = Ef1 + Ef2. With FULL model
%    this is true but with FIC and PIC this is true only
%    approximately. That is Ef \approx Ef1 + Ef2.
%
%    With CS+FIC the predictions are exact if the PREDCF covariance
%    functions are all in the FIC part or if they are CS
%    covariances.
%
%    NOTE! When making predictions with a subset of covariance
%    functions with FIC approximation the predictive variance can
%    in some cases be ill-behaved i.e. negative or
%    unrealistically small. This may happen because of the
%    approximative nature of the prediction.
%
%  See also
%    GP_SET, GP_OPTIM, DEMO_REGRESSION*
%

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2008 Jouni Hartikainen
% Copyright (c) 2011-2012 Ville Tolvanen
% Copyright (c) 2010,2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

if iscell(gp) || numel(gp.jitterSigma2)>1 || isfield(gp,'latent_method')
  % use an inference specific method
  if iscell(gp)
    fh_pred=@gpia_jpred;
  elseif numel(gp.jitterSigma2)>1
    fh_pred=@gpmc_jpred;
  elseif isfield(gp,'latent_method')
    fh_pred=gp.fh.jpred;
  else
    error('Logical error by coder of this function!')
  end
  switch nargout
    case 1
      [Eft] = fh_pred(gp, x, y, varargin{:});
    case 2
      [Eft, Covft] = fh_pred(gp, x, y, varargin{:});
    case 3
      [Eft, Covft, ljpyt] = fh_pred(gp, x, y, varargin{:});
    case 4
      [Eft, Covft, ljpyt, Eyt] = fh_pred(gp, x, y, varargin{:});
    case 5
      [Eft, Covft, ljpyt, Eyt, Covyt] = fh_pred(gp, x, y, varargin{:});
  end
  return
end

ip=inputParser;
ip.FunctionName = 'GP_JPRED';
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addOptional('xt', [], @(x) isempty(x) || (isreal(x) && all(isfinite(x(:)))))
ip.addParamValue('yt', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('predcf', [], @(x) isempty(x) || ...
                 isvector(x) && isreal(x) && all(isfinite(x)&x>0))
ip.addParamValue('tstind', [], @(x) isempty(x) || iscell(x) ||...
                 (isvector(x) && isreal(x) && all(isfinite(x)&x>0)))
ip.parse(gp, x, y, varargin{:});
xt=ip.Results.xt;
yt=ip.Results.yt;
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
end

tn = size(x,1);
if isfield(gp.lik, 'nondiagW') && ~ismember(gp.lik.type, {'LGP' 'LGPC'})
  switch gp.lik.type
    case {'Softmax', 'Multinom'}
      nout=size(y(:),1)/tn;
    otherwise
      if ~isfield(gp.lik,'xtime') && size(y,1)~=size(x,1)
        y=reshape(y,size(x,1),size(y,1)/size(x,1));
      end
      nout=length(gp.comp_cf);
      
      % Indices for looping over latent processes
      if ~isfield(gp.lik, 'xtime')
        nl=[0 repmat(n,1,nout)];
        nl=cumsum(nl);
      else
        xtime=gp.lik.xtime;
        ntime = size(xtime,1);
        n=tn-ntime;
        nl=[0 ntime n];
        nl=cumsum(nl);
      end
  end
  y=reshape(y,tn,nout);
  
  if isfield(gp, 'comp_cf')  % own covariance for each ouput component
    multicf = true;
    if length(gp.comp_cf) ~= nout
      error('GP2_PRED: the number of component vectors in gp.comp_cf must be the same as number of outputs or latent processes.')
    end
    if ~isempty(predcf)
      if ~iscell(predcf) || length(predcf)~=nout
        error(['GP2_PRED: if own covariance for each output or latent process component is used,'...
          'predcf has to be cell array and contain nout (vector) elements.   '])
      end
    else
      predcf = gp.comp_cf;
    end
  else
    multicf = false;
    for i1=1:nout
      predcf2{i1} = predcf;
    end
    predcf=predcf2;
  end
end

if nargout > 2 && isempty(yt)
  ljpyt=[];
end

% Evaluate this if sparse model is used
switch gp.type
  case 'FULL'
      
    if ~isfield(gp.lik, 'nondiagW') && ismember(gp.lik.type, {'LGP', 'LGPC'})
        %evaluate a = C\y;
        % -------------------
        [c, C]=gp_trcov(gp,x);
        
        if issparse(C)
          LD = ldlchol(C);
          a = ldlsolve(LD,y);
        elseif isempty(C)
          C=0;
          L=[];
          a = zeros(length(y),1);
        else
          L = chol(C)';
          a = L'\(L\y);
        end
        
        % evaluate K*a
        % -------------------
        K=gp_cov(gp,x,xt,predcf);
        Eft = K'*a;
        
        if  isfield(gp,'meanf')
          if issparse(C)
            [RB RAR] = mean_jpredf(gp,x,xt,K,LD,a,'gaussian',[]);    % terms with non-zero mean -prior
          else
            [RB RAR] = mean_jpredf(gp,x,xt,K,L,a,'gaussian',[]);    % terms with non-zero mean -prior
          end
          Eft = Eft + RB;
        end
        
        % Evaluate variance
        % Vector of diagonal elements of covariance matrix
        if nargout > 1
          
          V = gp_trcov(gp,xt,predcf);
          if issparse(C)
            Covft = (V - (K'*ldlsolve(LD,K)));
          else
            v = L\K;
            Covft = (V - (v'*v));
          end
          
          % If there are specified mean functions
          if  isfield(gp,'meanf')
            Covft = Covft + RAR;
          end
        end
        
        if nargout > 2
          % Scale mixture model in lik_smt is a special case
          % handle it separately
          if ~strcmp(gp.lik.type, 'lik_smt')
            % normal case
            [V, Cv] = gp_trvar(gp,xt,predcf);
            Eyt = Eft;
            apu = Cv - V;
            Covyt = Covft + diag(apu); % Utilize the Covft calculated above (faster!?) dimensions did not match here earlier!
            % Log joint predictive density
            if ~isempty(yt)
              ljpyt = mnorm_lpdf(yt', Eyt', Covyt);
            end
          else
            % scale mixture case
            nu = gp.lik.nu;             % Not working at the moment. Probably.
            sigma2 = gp.lik.tau2.*gp.lik.alpha.^2;
            sigma = sqrt(sigma2);
            
            Eyt = Eft;
            Covyt = (nu./(nu-2).*sigma2);
            
            for i2 = 1:length(Eft)
              mean_app = Eft(i2);
              sigm_app = sqrt(Covft(i2));
              
              pd = @(f) t_pdf(yt(i2), nu, f, sigma).*norm_pdf(f,Eft(i2),sqrt(Covft(i2)));
              pyt(i2) = quadgk(pd, mean_app - 12*sigm_app, mean_app + 12*sigm_app);
            end
          end
        end
    else
      % Likelihoods with non-diagonal Hessian
      L = zeros(tn,tn,nout);
      ntest=size(xt,1);
      K_nf = zeros(ntest,tn,nout);
      if multicf
        for i1=1:nout
          [tmp,C] = gp_trcov(gp, x, gp.comp_cf{i1});
          L(:,:,i1) = chol(C)';
          K_nf(:,:,i1) = gp_cov(gp,xt,x,predcf{i1});
        end
      else
        for i1=1:nout
          [tmp,C] = gp_trcov(gp, x);
          L(:,:,i1) = chol(C)';
          K_nf(:,:,i1) = gp_cov(gp,xt,x,predcf{i1});
        end
      end
            
      Eft = zeros(ntest,nout);
      for i1=1:nout
        Eft(:,i1) = K_nf(:,:,i1)*(L(:,:,i1)'\(L(:,:,i1)\y(:,i1)));
      end
      
      Covft = zeros(ntest*nout,ntest*nout);
      if nargout > 1
        if multicf
          for i1=1:nout
            v = L(:,:,i1)\K_nf(:,:,i1)';
            V = gp_trcov(gp,xt,gp.comp_cf{i1});
            Covft((i1-1)*ntest+1:i1*ntest,(i1-1)*ntest+1:i1*ntest) = V - v'*v;
          end
        else
          for i1=1:nout
            v = L(:,:,i1)\K_nf(:,:,i1)';
            V = gp_trcov(gp,xt,predcf{i1});
            Covft((i1-1)*ntest+1:i1*ntest,(i1-1)*ntest+1:i1*ntest) = V - v'*v;
          end
        end          
      end
      if nargout > 2
        % normal case
        Eyt = Eft;
        Covyt=Covft;
        ljpyt=0;
        for i1=1:nout
          [V, Cv] = gp_trvar(gp,xt,predcf{i1});
          % Diagonal indices
          i2=(i1-1)*(nout*ntest^2+ntest)+1:nout*ntest+1:i1*(nout*ntest^2+ntest);
          Covyt(i2) = Covyt(i2) + Cv' - V';
          if ~isempty(yt)
            ljpyt = ljpyt + mnorm_lpdf(yt(:,i1)', Eyt(:,i1)', Covyt((i1-1)*ntest+1:i1*ntest,(i1-1)*ntest+1:i1*ntest));
          end
        end
        Eyt=Eyt(:);
      end
      Eft=Eft(:);
    end

  case 'FIC'
    % Check the tstind vector
    if nargin > 5
        if ~isempty(tstind) && length(tstind) ~= size(x,1)
            error('tstind (if provided) has to be of same lenght as x.')
        end
    else
        tstind = [];
    end
    
    u = gp.X_u;
    m = size(u,1);
    % Turn the inducing vector on right direction
    if size(u,2) ~= size(x,2)
        u=u';
    end
    % Calculate some help matrices
    [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
    K_fu = gp_cov(gp, x, u);   % f x u
    K_uu = gp_trcov(gp, u);     % u x u, noiseles covariance K_uu
    K_nu = gp_cov(gp,xt,u);       % n x u
    Luu = chol(K_uu,'lower');
    
    % Evaluate the Lambda (La) for specific model
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    % Here we need only the diag(Q_ff), which is evaluated below
    B=Luu\(K_fu');
    Qv_ff=sum(B.^2)';
    Lav = Cv_ff-Qv_ff;   % 1 x f, Vector of diagonal elements
                         % iLaKfu = diag(inv(Lav))*K_fu = inv(La)*K_fu
    iLaKfu = zeros(size(K_fu));  % f x u,
    n=size(x,1)
    for i=1:n
        iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
    end
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;

    L = iLaKfu/chol(A);
    p = y./Lav - L*(L'*y);

    % Prediction matrices formed with only subset of cf's.
    if ~isempty(predcf)
        K_fu = gp_cov(gp, x, u, predcf);   % f x u
        K_uu = gp_trcov(gp, u, predcf);     % u x u, noiseles covariance K_uu
        K_nu = gp_cov(gp,xt,u,predcf);       % n x u
    end
    Eft = K_nu*(K_uu\(K_fu'*p));

    % if the prediction is made for training set, evaluate Lav also for prediction points
    if ~isempty(tstind)
        [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf);
        Luu = chol(K_uu,'lower');
        B=Luu\(K_fu');
        Qv_ff=sum(B.^2)';
        Lav2 = zeros(size(Eft));
        Lav2(tstind) = Kv_ff-Qv_ff;
        Eft(tstind) = Eft(tstind) + Lav2(tstind).*p;
    end

    if nargout > 1
        [Knn_v, Cnn_v] = gp_trvar(gp,xt,predcf);             
        Luu = chol(K_uu, 'lower');
        B=Luu\(K_fu');
        B2=Luu\(K_nu');
        Lav_n = Knn_v - sum(B2.^2)';
        BL = B*L;

        
        if isempty(tstind)
            error('tstind not provided.')
        end

        
        % if the prediction is made for training set, evaluate Lav also for prediction points
        if ~isempty(tstind)
          K = B'*B2 + diag(Lav_n);
          K2 = B2'*B2 + diag(Lav_n);
          C = B'*B + diag(Lav);
          
          L = chol(C,'lower');
          a = L'\(L\y);
          v = L\K;
          
          Covft = K2-v'*v;                         
        end
    end
        
    
    if nargout > 2
        Eyt = Eft;
        Covyt = Covft + diag(Cnn_v) - diag(Knn_v);
        if ~isempty(yt)
          ljpyt = mnorm_lpdf(yt', Eyt', Covyt);
        end
    end
    
  case {'PIC' 'PIC_BLOCK'}
    u = gp.X_u;
    ind = gp.tr_index;
    if size(u,2) ~= size(x,2)
        u=u';
    end

    % Calculate some help matrices
    [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
    K_fu = gp_cov(gp, x, u);         % f x u
    K_nu = gp_cov(gp, xt, u);         % n x u
    K_uu = gp_trcov(gp, u);    % u x u, noiseles covariance K_uu
    Luu = chol(K_uu)';

    % Evaluate the Lambda (La) for specific model
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    % Here we need only the diag(Q_ff), which is evaluated below
    B=Luu\K_fu';
    iLaKfu = zeros(size(K_fu));  % f x u
    for i=1:length(ind)
        Qbl_ff = B(:,ind{i})'*B(:,ind{i});
        [Kbl_ff, Cbl_ff] = gp_trcov(gp, x(ind{i},:));
        La{i} = Cbl_ff - Qbl_ff;
        iLaKfu(ind{i},:) = La{i}\K_fu(ind{i},:);    
    end
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;            % Ensure symmetry
    L = iLaKfu/chol(A);

    tyy = y;
    % From this on evaluate the prediction
    % See Snelson and Ghahramani (2007) for details
    p=iLaKfu*(A\(iLaKfu'*tyy));
    for i=1:length(ind)
        p2(ind{i},:) = La{i}\tyy(ind{i},:);
    end
    p= p2-p;
    
    % Prediction matrices formed with only subsetof cf's.
    if ~isempty(predcf)
        K_fu = gp_cov(gp, x, u, predcf);        % f x u
        K_nu = gp_cov(gp, xt, u, predcf);         % n x u
        K_uu = gp_trcov(gp, u, predcf);          % u x u, noiseles covariance K_uu
    end
        
    iKuuKuf = K_uu\K_fu';    
    w_bu=zeros(length(xt),length(u));
    w_n=zeros(length(xt),1);
    B2=Luu\(K_nu');
    for i=1:length(ind)
        w_bu(tstind{i},:) = repmat((iKuuKuf(:,ind{i})*p(ind{i},:))', length(tstind{i}),1);
        K_nf = gp_cov(gp, xt(tstind{i},:), x(ind{i},:),predcf);              % n x u
        w_n(tstind{i},:) = K_nf*p(ind{i},:);
        Qbl_ff = B2(:,tstind{i})'*B2(:,tstind{i});
        [Kbl_ff, Cbl_ff] = gp_trcov(gp, xt(tstind{i},:));
        La2{i} = Kbl_ff - Qbl_ff;
    end
    
    Eft = K_nu*(iKuuKuf*p) - sum(K_nu.*w_bu,2) + w_n;
    

    if nargout > 1        
        % Form iLaKfu again if a subset of cf's is used for making predictions
        if ~isempty(predcf)
            iLaKfu = zeros(size(K_fu));  % f x u
            for i=1:length(ind)
                iLaKfu(ind{i},:) = La{i}\K_fu(ind{i},:);    
            end
        end
        
%         kstarstar = gp_trvar(gp, xt, predcf);       
%         KnuiKuu = K_nu/K_uu;
%         KufiLaKfu = K_fu'*iLaKfu;
%         QnfL = KnuiKuu*(K_fu'*L);
%         Covft1 = zeros(size(xt,1),size(xt,1));
%         Covft2 = zeros(size(xt,1),size(xt,1));
%         Covft3 = zeros(size(xt,1),size(xt,1));
%         C = B'*B;
%         KKnn = gp_trcov(gp,xt,predcf);
%         Knn = B2'*B2;
%         Knf = B2'*B;
%         KKnf = gp_cov(gp,xt,x,predcf);
%         for i=1:length(ind)
%             KubiLaKbu = K_fu(ind{i},:)'/La{i}*K_fu(ind{i},:);
%             nonblock = KufiLaKfu - KubiLaKbu;
%             Covft1(tstind{i}) = diag(KnuiKuu(tstind{i},:)*nonblock*KnuiKuu(tstind{i},:)');
%             
%             Knb = gp_cov(gp, xt(tstind{i},:), x(ind{i},:), predcf);
%             Covft2(tstind{i}) = diag(Knb/La{i}*Knb');
%             
%             KnbL = Knb*L(ind{i},:);
%             QnbL = KnuiKuu(tstind{i},:)*(K_fu(ind{i},:)'*L(ind{i},:));
%             %Covft3(tstind{i}) = sum(QnfL(tstind{i},:) - QnbL + KnbL,2);
%             Covft3(tstind{i}) = diag((QnfL(tstind{i},:) - QnbL + KnbL)*(QnfL(tstind{i},:) - QnbL + KnbL)');
%             C(ind{i},ind{i}) = C(ind{i},ind{i}) + La{i};
%             Knn(ind{i},ind{i}) = Knn(ind{i},ind{i}) + KKnn(tstind{i},tstind{i}) - B2(:,tstind{i})'*B2(:,tstind{i});
%             Knf(tstind{i},ind{i}) = Knf(tstind{i},ind{i}) + KKnf(tstind{i},ind{i}) - B2(:,tstind{i})'*B(:,ind{i});
%         end
%         L = chol(C,'lower');
%         v = L\Knf;
%         Covft = Knn-v'*v;
% %         Covft = kstarstar - (Covft1 + Covft2 - Covft3);       %MUUTETTU
%     end
    
 
        
        B2=Luu\(K_nu');
        C = B'*B;
        KKnn = gp_trcov(gp,xt,predcf);
        Knn = B2'*B2;
        Knf = B2'*B;
        KKnf = gp_cov(gp,xt,x,predcf);
        for i=1:length(ind)
            if (size(ind{i},1) ~= size(tstind{i},1))
                error('Testing and training block cell array vectors differ in size.');     % gp.tr_index and tstind
            end
            C(ind{i},ind{i}) = C(ind{i},ind{i}) + La{i};
            Knn(ind{i},ind{i}) = Knn(ind{i},ind{i}) + KKnn(tstind{i},tstind{i}) - B2(:,tstind{i})'*B2(:,tstind{i});
            Knf(tstind{i},ind{i}) = Knf(tstind{i},ind{i}) + KKnf(tstind{i},ind{i}) - B2(:,tstind{i})'*B(:,ind{i});
        end

        L = chol(C)';

        v = L\Knf';
        Covft = (Knn) - (v'*v);
    end
    
    if nargout > 2
        Eyt = Eft;
        [Knn_v, Cnn_v] = gp_trvar(gp,xt,predcf);
        Covyt = Covft + diag(Cnn_v) - diag(Knn_v);
        if ~isempty(yt)
          ljpyt = mnorm_lpdf(yt', Eyt', Covyt);
        end
    end
  case 'CS+FIC'
    % Here tstind = 1 if the prediction is made for the training set 
    if nargin > 5
        if ~isempty(tstind) && length(tstind) ~= size(x,1)
            error('tstind (if provided) has to be of same lenght as x.')
        end
    else
        tstind = [];
    end
    
    n = size(x,1);
    n2 = size(xt,1);

    u = gp.X_u;
    m = size(u,1);
    ncf = length(gp.cf);
    
    % Indexes to all non-compact support and compact support covariances.
    cf1 = [];
    cf2 = [];
    % Indexes to non-CS and CS covariances, which are used for predictions
    predcf1 = [];
    predcf2 = [];    

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
    
    % First evaluate needed covariance matrices
    % v defines that parameter is a vector
    [Kv_ff, Cv_ff] = gp_trvar(gp, x, cf1);  % f x 1  vector    
    K_fu = gp_cov(gp, x, u, cf1);         % f x u
    K_uu = gp_trcov(gp, u, cf1);    % u x u, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu

    Luu  = chol(K_uu)';
    K_nu = gp_cov(gp, xt, u, cf1);         % n x u

    % Evaluate the Lambda (La)
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    B=Luu\(K_fu');       % u x f
    Qv_ff=sum(B.^2)';
    Lav = Cv_ff-Qv_ff;   % f x 1, Vector of diagonal elements

    K_cs = gp_trcov(gp,x,cf2);
    Kcs_nf = gp_cov(gp, xt, x, predcf2);
    Kcs_nn = gp_trcov(gp, xt, predcf2);
    La = sparse(1:tn,1:tn,Lav,tn,tn) + K_cs;
    
    iLaKfu = La\K_fu;
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;     % Ensure symmetry
    L = iLaKfu/chol(A);
    
    p = La\y - L*(L'*y);

    %p2 = y./Lav - iLaKfu*(A\(iLaKfu'*y));
    %    Knf = K_nu*(K_uu\K_fu');

    K_fu = gp_cov(gp, x, u, predcf1);       % f x u
    K_uu = gp_trcov(gp, u, predcf1);         % u x u, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;                  % ensure the symmetry of K_uu
    K_nu = gp_cov(gp, xt, u, predcf1);        % n x u    

    % Calculate the predictive mean according to the type of
    % covariance functions used for making the prediction
    if ptype == 1
        Eft = K_nu*(K_uu\(K_fu'*p));
    elseif ptype == 2
        Eft = Kcs_nf*p;
    else 
        Eft = K_nu*(K_uu\(K_fu'*p)) + Kcs_nf*p;
    end
    
    % evaluate also Lav2 if the prediction is made for training set
    if ~isempty(tstind)
        [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf1);
        Luu = chol(K_uu)';
        B=Luu\(K_fu');
        Qv_ff=sum(B.^2)';
        Lav2 = zeros(size(Eft));
        Lav2(tstind) = Kv_ff-Qv_ff;
    end  

    % Add also Lav2 if the prediction is made for training set
    % and non-CS covariance function is used for prediction
    if ~isempty(tstind) && (ptype == 1 || ptype == 3)
        Eft(tstind) = Eft(tstind) + Lav2(tstind).*p;
    end
    
    if nargout > 1
        [Knn_v, Cnn_v] = gp_trvar(gp,xt,predcf);
        Luu = chol(K_uu, 'lower');
        B=Luu\(K_fu');
        B2=Luu\(K_nu');
        iLaKfu = La\K_fu;
        
        % Calculate the predictive variance according to the type
        % covariance functions used for making the prediction
        if ptype == 1 || ptype == 3                            
            % FIC part of the covariance
%             Covft = Knn_v - sum(B2'.*(B*(La\B')*B2)',2) +
%             sum((K_nu*(K_uu\(K_fu'*L))).^2, 2);             % MUUTETTU
            tmpmatr = (K_nu*(K_uu\(K_fu'*L)));
            Covft = B2'*B2 + Kcs_nn + sparse(1:n2,1:n2,Knn_v - sum(B2.^2)',n2,n2) - B2'*(B*(La\B')*B2) + tmpmatr*tmpmatr';
%             Covft = B2'*B2 + Kcs_nn - B2'*(B*(La\B')*B2) + tmpmatr*tmpmatr';


            % Add Lav2 if the prediction is made for the training set
            if  ~isempty(tstind)
                % Non-CS covariance
                if ptype == 1
                    Kcs_nf = sparse(tstind,1:n,Lav2(tstind),n2,n);
                % Non-CS and CS covariances
                else
                    Kcs_nf = Kcs_nf + sparse(tstind,1:n,Lav2(tstind),n2,n);
                end
                 % Add Lav2 and possibly Kcs_nf
                 tmpmatr = Kcs_nf/chol(La);
                 Covft = Covft - tmpmatr*tmpmatr';
                 tmpmatr = Kcs_nf*L;
                 Covft = Covft + tmpmatr*tmpmatr';
                 Covft = Covft - 2.*(Kcs_nf*iLaKfu)*(K_uu\K_nu') + 2.*(Kcs_nf*L)*(L'*K_fu*(K_uu\K_nu'));
            % In case of both non-CS and CS prediction covariances add 
            % only Kcs_nf if the prediction is not done for the training set 
            elseif ptype == 3
                tmpmatr = Kcs_nf/chol(La);
                Covft = Covft - tmpmatr*tmpmatr';
                tmpmatr = Kcs_nf*L;
                Covft = Covft + tmpmatr*tmpmatr';
                Covft = Covft - 2*(Kcs_nf*iLaKfu)*(K_uu\K_nu') + 2.*(Kcs_nf*L)*(L'*K_fu*(K_uu\K_nu'));
            end
        % Prediction with only CS covariance
        elseif ptype == 2
            Covft = Knn_v - sum((Kcs_nf/chol(La)).^2,2) + sum((Kcs_nf*L).^2, 2) ;

        end        
    end
    
% $$$     Lav_pr = Kv_ff-Qv_ff;
% $$$     K2 = B'*B + K_cs + diag(Lav_pr);
% $$$     C = B'*B + K_cs + diag(Lav);
% $$$     K = B'*B + Kcs_nf + diag(Lav_pr); 
% $$$     
% $$$     L = chol(C)';
% $$$     %    y=K'*(C\y);
% $$$     a = L'\(L\y);
% $$$     Eft = K'*a;
% $$$ 
% $$$     v = L\K;
% $$$     V = gp_trvar(gp,xt,predcf);
% $$$     Covft = V - diag(v'*v);

    if nargout > 2
        Eyt = Eft;
        Covyt = Covft + diag(Cnn_v) - diag(Knn_v);
        if ~isempty(yt)
          ljpyt = mnorm_lpdf(yt', Eyt', Covyt);
        end
    end
    
  case {'VAR' 'DTC' 'SOR'}
    % Check the tstind vector
    if nargin > 5
        if ~isempty(tstind) && length(tstind) ~= size(tx,1)
            error('tstind (if provided) has to be of same lenght as tx.')
        end
    else
        tstind = [];
    end
    
    u = gp.X_u;
    m = size(u,1);
    % Turn the inducing vector on right direction
    if size(u,2) ~= size(x,2)
        u=u';
    end
    % Calculate some help matrices
    [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
    K_fu = gp_cov(gp, x, u);   % f x u
    K_uu = gp_trcov(gp, u);     % u x u, noiseles covariance K_uu
    K_nu = gp_cov(gp,xt,u);       % n x u
    Luu = chol(K_uu, 'lower');
    
    % Evaluate the Lambda (La) for specific model
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    % Here we need only the diag(Q_ff), which is evaluated below
    B=Luu\(K_fu');
    Qv_ff=sum(B.^2)';
    Lav2 = diag(Cv_ff-Kv_ff);
    Lav = Cv_ff-Kv_ff;   % 1 x f, Vector of diagonal elements
                         % iLaKfu = diag(inv(Lav))*K_fu = inv(La)*K_fu
    iLaKfu = zeros(size(K_fu));  % f x u,
    n=size(x,1)
    for i=1:n
        iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
    end
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;

    L = iLaKfu/chol(A);
    p = y./Lav - L*(L'*y);

    % Prediction matrices formed with only subset of cf's.
    if ~isempty(predcf)
        K_fu = gp_cov(gp, x, u, predcf);   % f x u
        K_uu = gp_trcov(gp, u, predcf);     % u x u, noiseles covariance K_uu
        K_nu = gp_cov(gp,xt,u,predcf);       % n x u
    end
    Eft = K_nu*(K_uu\(K_fu'*p));


    if nargout > 1
%         [Knn_v, Cnn_v] = gp_trvar(gp,xt,predcf);
        [Knn_v, Cnn_v] = gp_trcov(gp,xt,predcf);
        Luu = chol(K_uu,'lower');
        B=Luu\(K_fu');
        B2=Luu\(K_nu');
        
        Covftr = sum(B2'.*(B*bsxfun(@ldivide,Lav2,B')*B2)',2) - sum((K_nu*(K_uu\(K_fu'*L))).^2, 2);
%         Covftr = 
        switch gp.type
          case {'VAR' 'DTC'}
            Covft = Knn_v - Covftr;         % Knn_v = diag(K_*,*)
          case  'SOR'
            Covft = sum(B2.^2,1)' - Covftr;     % sum(B2.^2,1)' = diag(Q_*,*)
        end

    end
    if nargout > 2
        Eyt = Eft;
        switch gp.type
          case {'VAR' 'DTC'}
            Covyt = Covft + Cnn_v - Knn_v;
          case 'SOR'
            Covyt = Covft + Cnn_v - sum(B2.^2,1)';
        end
    end
    if nargout > 4
        pyt = norm_pdf(y, Eyt, sqrt(Covyt));
    end  
    
  case 'SSGP'
    if nargin > 4
        error(['Prediction with a subset of original ' ...
               'covariance functions not currently implemented with SSGP']);
    end

    [Phi_f, S] = gp_trcov(gp, x);
    Phi_a = gp_trcov(gp, xt);
    m = size(Phi_f,2);
    ns = eye(m,m)*S(1,1);
    
    L = chol(Phi_f'*Phi_f + ns,'lower');
    Eft = Phi_a*(L'\(L\(Phi_f'*y)));

    
    if nargout > 1
        Covft = sum(Phi_a/L',2)*S(1,1);
    end
    if nargout > 2
        error('gp_pred with three output arguments is not implemented for SSGP!')
    end
end