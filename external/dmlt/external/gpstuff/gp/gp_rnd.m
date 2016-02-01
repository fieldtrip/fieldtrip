function [sampft, sampyt] = gp_rnd(gp, x, y, xt, varargin)
%GP_RND  Random draws from the posterior Gaussian process
%
%  Description
%    [SAMPFT, SAMPYT] = GP_RND(GP, X, Y, XT, OPTIONS) takes a
%    Gaussian process structure, record structure (from gp_mc) or
%    array (from gp_ia) GP together with a matrix XT of input
%    vectors, matrix X of training inputs and vector Y of training
%    targets, and returns a random sample SAMPFT and SAMPYT from
%    the posterior distribution p(ft|x,y,xt) and the predictive
%    distribution p(yt|x,y,xt) at locations XT.
%
%    OPTIONS is optional parameter-value pair
%      nsamp  - determines the number of samples (default = 1).
%      predcf - index vector telling which covariance functions are 
%               used for prediction. Default is all (1:gpcfn)
%      tstind - a vector defining, which rows of X belong to which 
%               training block in *IC type sparse models. Default is [].
%               See also GP_PRED.
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case. 
%
%    If likelihood is non-Gaussian and gp.latent_method is either
%    Laplace or EP, the samples are drawn from the Gaussian
%    posterior approximation obtained from gpla_e or gpep_e.
%
%  See also
%    GP_PRED, GP_PAK, GP_UNPAK

%  Internal comments
%    - sampling with FIC, PIC and CS+FIC forms full nxn matrix and
%       works only when sampling for the training inputs
%    - The code is not optimized

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2008      Jouni Hartikainen
% Copyright (c) 2011      Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
ip.FunctionName = 'GP_RND';
ip.addRequired('gp',@(x) isstruct(x) || iscell(x));
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('xt', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('zt', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('predcf', [], @(x) isempty(x) || ...
                 isvector(x) && isreal(x) && all(isfinite(x)&x>0))
ip.addParamValue('tstind', [], @(x) isempty(x) || iscell(x) ||...
                 (isvector(x) && isreal(x) && all(isfinite(x)&x>0)))
ip.addParamValue('nsamp', 1, @(x) isreal(x) && isscalar(x))
ip.parse(gp, x, y, xt, varargin{:});
z=ip.Results.z;
zt=ip.Results.zt;
predcf=ip.Results.predcf;
tstind=ip.Results.tstind;
nsamp=ip.Results.nsamp;

tn = size(x,1);

if isstruct(gp) && numel(gp.jitterSigma2)==1
  % Single GP
  if isfield(gp.lik.fh,'trcov') || isfield(gp, 'latentValues')
    % ===================================
    % Gaussian likelihood or MCMC with latent values
    % ===================================
    
    % Evaluate this if sparse model is used
    switch gp.type
      case 'FULL'
        [c, C]=gp_trcov(gp,x);
        K=gp_cov(gp,x,xt,predcf);
        [K2, C2] = gp_trcov(gp,xt,predcf);
                
        if issparse(C)
          LD = ldlchol(C);
          Ef = repmat( K'*ldlsolve(LD,y), 1, nsamp) ;
          pcov = K2 - K'*ldlsolve(LD,K);
          if  isfield(gp,'meanf')
              [RB RAR] = mean_jpredf(gp,x,xt,K,LD,a,'gaussian',[]);    % terms with non-zero mean -prior
              Ef = Ef + repmat(RB,1,nsamp);
              pcov = pcov + RAR;
          end
          predcov = chol(pcov,'lower');
          sampft = Ef + predcov*randn(size(Ef));
          if nargout > 1
            pcov = C2 - K'*ldlsolve(LD,K);
            if  isfield(gp,'meanf')
                pcov = pcov + RAR;
            end
            predcov = chol(pcov,'lower');
            sampyt = Ef + predcov*randn(size(Ef));
          end        
        else
          L = chol(C,'lower');
          %    y=K'*(C\y);
          a = L'\(L\y);
          Ef = repmat( K'*a, 1, nsamp);
          v = L\K;

          pcov = K2-v'*v;
          if  isfield(gp,'meanf')
              [RB RAR] = mean_jpredf(gp,x,xt,K,L,a,'gaussian',[]);    % terms with non-zero mean -prior
              Ef = Ef + repmat(RB,1,nsamp);
              pcov = pcov + RAR;
          end
          predcov = chol(pcov,'lower');
          sampft = Ef + predcov*randn(size(Ef));
          if nargout > 1
            pcov = C2-v'*v;
            if  isfield(gp,'meanf')
                pcov = pcov + RAR;
            end
            predcov = chol(pcov,'lower');
            sampyt = Ef + predcov*randn(size(Ef));
          end
        end   
        
      case 'FIC'    
        % Here tstind = 1 if the prediction is made for the training set 
        if nargin > 5
          if length(tstind) ~= size(x,1) && ~isempty(tstind)
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
        n=size(x,1);
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
        Ef = K_nu*(K_uu\(K_fu'*p)) ;

        % if the prediction is made for training set, evaluate Lav also for
        % prediction points
        if ~isempty(tstind)
          [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf);
          Luu = chol(K_uu,'lower');
          B=Luu\(K_fu');
          Qv_ff=sum(B.^2)';
          Lav2 = zeros(size(Ef));
          Lav2(tstind) = Kv_ff-Qv_ff;
          Ef(tstind) = Ef(tstind) + Lav2(tstind).*p;
        end
        
        % Sigma_post = Qnn + La_n - Qnf*(Qff+La_f)^(-1)*Qfn
        %            = B2'*(I-B*La_f^(-1)*B' + B*L*L'*B')*B2 + La_n 
        %
        % in case of tstind is given:
        % Sigma_post = Qnn + La_n - (Qnf+la2)*(Qff+La_f)^(-1)*(Qfn+La2)

        
        Ef = repmat(Ef , 1, nsamp);
        
        [Knn_v, Cnn_v] = gp_trvar(gp,xt,predcf);
        Luu = chol(K_uu,'lower');
        B=Luu\(K_fu');
        B2 = Luu\(K_nu');
        Lav_n = Knn_v - sum(B2.^2)';
        BL = B*L;

        Sigm_mm = eye(size(K_uu)) - B*(repmat(Lav,1,m).\B') + BL*BL';
        sampft = Ef + B2'*(chol(Sigm_mm)'*randn(m,nsamp)) + randn(size(Ef)).*sqrt(repmat(Lav_n,1,nsamp));

        if ~isempty(tstind)
          K = B'*B2 + diag(Lav_n);
          K2 = B2'*B2 + diag(Lav_n);
          C = B'*B + diag(Lav);
          
          L = chol(C,'lower');
          %    y=K'*(C\y);
          a = L'\(L\y);
          Ef = repmat( K'*a, 1, nsamp);
          v = L\K;
          
          predcov = chol(K2-v'*v,'lower');
          sampft = Ef + predcov*randn(size(Ef));        
        end

        if nargout > 1
          sigma = sqrt(Cnn_v-Knn_v);
          sampyt = sampft + sigma*randn(size(sampft));
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
        p = p2-p;
        
        % Prediction matrices formed with only subsetof cf's.
        if ~isempty(predcf)
          K_fu = gp_cov(gp, x, u, predcf);        % f x u
          K_nu = gp_cov(gp, xt, u, predcf);         % n x u
          K_uu = gp_trcov(gp, u, predcf);          % u x u, noiseles covariance K_uu
        end
        
        iKuuKuf = K_uu\K_fu';    
        w_bu=zeros(length(xt),length(u));
        w_n=zeros(length(xt),1);
        B2 = Luu\(K_nu');
        for i=1:length(ind)
          w_bu(tstind{i},:) = repmat((iKuuKuf(:,ind{i})*p(ind{i},:))', length(tstind{i}),1);
          K_nf = gp_cov(gp, xt(tstind{i},:), x(ind{i},:),predcf);              % n x u
          w_n(tstind{i},:) = K_nf*p(ind{i},:);
          
          Qbl_ff = B2(:,tstind{i})'*B2(:,tstind{i});
          [Kbl_ff, Cbl_ff] = gp_trcov(gp, xt(tstind{i},:));
          La2{i} = Kbl_ff - Qbl_ff;
          La22{i} = Cbl_ff - Qbl_ff;
        end
        
        Ef = repmat(K_nu*(iKuuKuf*p) - sum(K_nu.*w_bu,2) + w_n, 1, nsamp);
        
        % Sigma_post = Qnn + La_n - Qnf*(Qff+La_f)^(-1)*Qfn
        %            = B'*(I-B*La_f^(-1)*B' + B*L*L'*B')*B + La_n
        BL = B*L;
        sampft = randn(size(Ef));
        for i=1:length(ind)
          iLaB(ind{i},:) = La{i}\B(:,ind{i})';
          sampft(ind{i},:) = chol(La2{i})'*sampft(ind{i},:);
        end
        Sigm_mm = eye(size(K_uu)) - B*iLaB + BL*BL';
        
        sampft = Ef + B2'*(chol(Sigm_mm)'*randn(size(K_uu,1),nsamp)) + sampft;
        
        if ~isempty(tstind)
          K = B'*B2;
          K2 = B2'*B2;
          C = B'*B;
          
          for i = 1:length(ind)
            K(ind{i},tstind{i}) = K(ind{i},tstind{i}) + La2{i};
            K2(tstind{i},tstind{i}) = K2(tstind{i},tstind{i}) + La2{i};
            C(tstind{i},tstind{i}) = C(ind{i},ind{i}) + La{i};
          end
          L = chol(C)';
          %    y=K'*(C\y);
          a = L'\(L\y);
          Ef = repmat( K'*a, 1, nsamp);
          v = L\K;
          
          predcov = chol(K2-v'*v, 'lower');
          sampft = Ef + predcov*randn(size(Ef));        
        end

        if nargout > 1
          sigma = sqrt(Cnn_v-Knn_v);
          sampyt = sampft + sigma*randn(size(sampft));
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
          Ef = K_nu*(K_uu\(K_fu'*p));
        elseif ptype == 2
          Ef = Kcs_nf*p;
        else 
          Ef = K_nu*(K_uu\(K_fu'*p)) + Kcs_nf*p;        
        end
        
        % evaluate also Lav2 if the prediction is made for training set
        if ~isempty(tstind)
          [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf1);
          Luu = chol(K_uu)';
          B=Luu\(K_fu');
          Qv_ff=sum(B.^2)';
          Lav2 = zeros(size(Ef));
          Lav2(tstind) = Kv_ff-Qv_ff;
        end  

        % Add also Lav2 if the prediction is made for training set
        % and non-CS covariance function is used for prediction
        if ~isempty(tstind) && (ptype == 1 || ptype == 3)
          Ef(tstind) = Ef(tstind) + Lav2(tstind).*p(tstind);
        end

        
        [Knn_v, Cnn_v] = gp_trvar(gp,xt,predcf1);
        Luu = chol(K_uu)';
        B=Luu\(K_fu');
        B2=Luu\(K_nu');
        iLaKfu = La\K_fu;
        
        % Calculate the predictive variance according to the type
        % covariance functions used for making the prediction
        
        % NOTE! We form full matrices below. Needs to be rewritten so that at most nxm matrices are formed
        
        if ptype == 1 || ptype == 3                            
          % FIC part of the covariance
          tmpmatr = (K_nu*(K_uu\(K_fu'*L)));
          Covf = B2'*B2 + Kcs_nn + sparse(1:n2,1:n2,Knn_v - sum(B2.^2)',n2,n2) - B2'*(B*(La\B')*B2) + tmpmatr*tmpmatr';
          
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
            Covf = Covf - tmpmatr*tmpmatr';
            tmpmatr = Kcs_nf*L;
            Covf = Covf + tmpmatr*tmpmatr';
            Covf = Covf - 2.*(Kcs_nf*iLaKfu)*(K_uu\K_nu') + 2.*(Kcs_nf*L)*(L'*K_fu*(K_uu\K_nu'));

            % In case of both non-CS and CS prediction covariances add 
            % only Kcs_nf if the prediction is not done for the training set 
          elseif ptype == 3
            tmpmatr = Kcs_nf/chol(La);
            Covf = Covf - tmpmatr*tmpmatr';
            tmpmatr = Kcs_nf*L;
            Covf = Covf + tmpmatr*tmpmatr';
            Covf = Covf - 2*(Kcs_nf*iLaKfu)*(K_uu\K_nu') + 2.*(Kcs_nf*L)*(L'*K_fu*(K_uu\K_nu'));
          end
          % Prediction with only CS covariance
        elseif ptype == 2
          Covf = Knn_v - sum((Kcs_nf/chol(La)).^2,2) + sum((Kcs_nf*L).^2, 2) ;
        end        
        
        Ef = repmat(Ef, 1, nsamp);
        predcov = chol(Covf,'lower');
        sampft = Ef + predcov*randn(size(Ef));
        
        if nargout > 1
          sigma = sqrt(Cnn_v-Knn_v);
          sampyt = sampft + sigma*randn(size(sampft));
        end
        
% $$$     % Check results with full matrices    
% $$$     Lav_pr = Kv_ff-Qv_ff;
% $$$     
% $$$     K2 = B2'*B2 + Kcs_nn + diag(Knn_v - sum(B2.^2)');
% $$$     C = B'*B + K_cs + diag(Lav);
% $$$     K = B2'*B + K_cs + diag(Lav_pr); 
% $$$     
% $$$     L = chol(C)';
% $$$     %    y=K'*(C\y);
% $$$     a = L'\(L\y);
% $$$     Ef = repmat( K'*a, 1, nsamp);
% $$$     v = L\K;
% $$$     
% $$$     predcov = chol(K2-v'*v)';
% $$$     sampft = Ef + predcov*randn(size(Ef));
        
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
        Ef = Phi_a*(L'\(L\(Phi_f'*y)));

        
        if nargout > 1
          Covf = sum(Phi_a/L',2)*S(1,1);
        end
        if nargout > 2
          error('gp_pred with three output arguments is not implemented for SSGP!')
        end
    end
    

  else
    % ===================================
    % non-Gaussian likelihood
    % ===================================
    switch gp.type
      % ---------------------------
      case 'FULL'

        switch gp.latent_method
          case 'Laplace'
            [e, edata, eprior, f, L] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
            
            W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
            deriv = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
            ntest=size(xt,1);
            
            % Evaluate the expectation
            K_nf = gp_cov(gp,xt,x,predcf);
            Ef = K_nf*deriv;
            
            % Evaluate the variance
            K = gp_trcov(gp,xt,predcf);
            if W >= 0
              if issparse(K_nf) && issparse(L)
                K = gp_trcov(gp, x);
                sqrtW = sparse(1:tn, 1:tn, sqrt(W), tn, tn);
                sqrtWKfn = sqrtW*K_nf';
                V = ldlsolve(L,sqrtWKfn);
                Covf = K - sqrtWKfn'*V;
              else
                sW = diag(sqrt(W));
                V = L\(sW*K_nf');
                Covf = K - V'*V;
              end
            else
              V = L*diag(W);
              R = diag(W) - V'*V;
              Covf = K - K_nf*(R*K_nf');
            end
          case 'EP'
            
            [e, edata, eprior, tautilde, nutilde, L] = gpep_e(gp_pak(gp), gp, x, y, 'z', z);
            
            [K, C]=gp_trcov(gp,x);
            K = gp_trcov(gp, xt, predcf);
            ntest=size(xt,1);
            K_nf=gp_cov(gp,xt,x,predcf);
            [n,nin] = size(x);
            
            if all(tautilde > 0) && ~isequal(gp.latent_opt.optim_method, 'robust-EP')
              sqrttautilde = sqrt(tautilde);
              Stildesqroot = sparse(1:n, 1:n, sqrttautilde, n, n);
              
              if issparse(L)
                z=Stildesqroot*ldlsolve(L,Stildesqroot*(C*nutilde));
              else
                z=Stildesqroot*(L'\(L\(Stildesqroot*(C*nutilde))));
              end
              Ef=K_nf*(nutilde-z);

              % Compute variance
              if issparse(L)
                V = ldlsolve(L, Stildesqroot*K_nf');
                Covf = K - K_nf*(Stildesqroot*V);
              else
                V = (L\Stildesqroot)*K_nf';
                Covf = K - V'*V;
              end
            else
              z=tautilde.*(L'*(L*nutilde));
              Ef=K_nf*(nutilde-z);
              
              S = diag(tautilde);
              V = K_nf*S*L';
              Covf = K - (K_nf*S)*K_nf' + V*V';
            end
        end
        
        predcov = chol(Covf,'lower');
        Ef = repmat(Ef,1,nsamp);
        sampft = Ef + predcov*randn(size(Ef));
        
        % ---------------------------
      case 'FIC'
        
        switch gp.latent_method
          case 'Laplace'
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
            K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu
            Luu = chol(K_uu)';
            
            m = size(u,1);
            
            [e, edata, eprior, f, L, a, La2] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);

            deriv = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
            ntest=size(xt,1);
            
            K_nu=gp_cov(gp,xt,u,predcf);
            Ef = K_nu*(Luu'\(Luu\(K_fu'*deriv)));
            
            % if the prediction is made for training set, evaluate Lav also for prediction points
            if ~isempty(tstind)
              [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf);
              K_fu = gp_cov(gp, x, u);
              B=Luu\(K_fu');
              Qv_ff=sum(B.^2)';
              %Lav = zeros(size(La));
              %Lav(tstind) = Kv_ff-Qv_ff;
              Lav = Kv_ff-Qv_ff;            
              Ef(tstind) = Ef(tstind) + Lav.*deriv;
            end
            
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

            
            % NOTE!
            % This is done with full matrices at the moment. 
            % Needs to be rewritten.
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
            Covf = K - V'*V;
            
            predcov = chol(Covf,'lower');
            Ef = repmat(Ef,1,nsamp);
            sampft = Ef + predcov*randn(size(Ef));
            
          case 'EP'
            [e, edata, eprior, tautilde, nutilde, L, La, b] = gpep_e(gp_pak(gp), gp, x, y, 'z', z);

            % Here tstind = 1 if the prediction is made for the training set 
            if nargin > 6
              if ~isempty(tstind) && length(tstind) ~= size(x,1)
                error('tstind (if provided) has to be of same lenght as x.')
              end
            else
              tstind = [];
            end
            
            u = gp.X_u;
            m = size(u,1);
            
            K_fu = gp_cov(gp,x,u,predcf);         % f x u
            K_nu=gp_cov(gp,xt,u,predcf);
            K_uu = gp_trcov(gp,u,predcf);          % u x u, noiseles covariance K_uu
            K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu

            kstarstar=gp_trvar(gp,xt,predcf);        
            
            % From this on evaluate the prediction
            % See Snelson and Ghahramani (2007) for details 
            %        p=iLaKfu*(A\(iLaKfu'*myytilde));
            p = b';
            
            ntest=size(xt,1);
            
            Ef = K_nu*(K_uu\(K_fu'*p));
            
            % if the prediction is made for training set, evaluate Lav also for prediction points
            if ~isempty(tstind)
              [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf);
              Luu = chol(K_uu)';
              B=Luu\(K_fu');
              Qv_ff=sum(B.^2)';
              Lav = Kv_ff-Qv_ff;
              Ef(tstind) = Ef(tstind) + Lav.*p;
            end
            
            % NOTE!
            % This is done with full matrices at the moment. 
            % Needs to be rewritten.            
            Luu = chol(K_uu)';
            B=Luu\(K_fu');   
            B2=Luu\(K_nu');   

            Knf = B2'*B;
            Knn = B2'*B2 + diag(kstarstar - sum(B2.^2)');
            if ~isempty(tstind)
              Knf(tstind,:) = Knf(tstind,:) + diag(kstarstar(tstind) - sum(B.^2)');
            end
            
            Covf = Knn - Knf * ( diag(1./La) - L*L' ) * Knf';
            
            predcov = chol(Covf,'lower');
            Ef = repmat(Ef,1,nsamp);
            sampft = Ef + predcov*randn(size(Ef));
            
        end        

        % ---------------------------
      case 'PIC'
        
        u = gp.X_u;
        K_fu = gp_cov(gp, x, u, predcf);         % f x u
        K_uu = gp_trcov(gp, u, predcf);          % u x u, noiseles covariance K_uu
        K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu
        K_nu=gp_cov(gp,xt,u,predcf);
        
        ind = gp.tr_index;
        ntest = size(xt,1);
        m = size(u,1);
        Luu = chol(K_uu)';
        B=Luu\(K_fu');
        B2 = Luu\(K_nu');

        switch gp.latent_method
          case 'Laplace'
            
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
            
            Ef = K_nu*(iKuuKuf*deriv) - sum(K_nu.*w_bu,2) + w_n;
            
            % Evaluate the variance
            W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
            kstarstar = gp_trvar(gp,xt,predcf);
            sqrtW = sqrt(W);
            % Components for (I + W^(1/2)*(Qff + La2)*W^(1/2))^(-1) = Lahat^(-1) - L2*L2'
            for i=1:length(ind)
              La{i} = diag(sqrtW(ind{i}))*La2{i}*diag(sqrtW(ind{i}));
              Lahat{i} = eye(size(La{i})) + La{i};
            end
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
            
            Covf = Knn - Knf * C * Knf';
            
            predcov = chol(Covf,'lower');
            Ef = repmat(Ef,1,nsamp);
            sampft = Ef + predcov*randn(size(Ef));

            
          case 'EP'
            
            [e, edata, eprior, tautilde, nutilde, L, La, b] = gpep_e(gp_pak(gp), gp, x, y, 'z', z);
            
            p = b';

            iKuuKuf = K_uu\K_fu';
            
            w_bu=zeros(length(xt),length(u));
            w_n=zeros(length(xt),1);
            for i=1:length(ind)
              w_bu(tstind{i},:) = repmat((iKuuKuf(:,ind{i})*p(ind{i},:))', length(tstind{i}),1);
              K_nf = gp_cov(gp, xt(tstind{i},:), x(ind{i},:), predcf);              % n x u
              w_n(tstind{i},:) = K_nf*p(ind{i},:);
            end
            
            Ef = K_nu*(iKuuKuf*p) - sum(K_nu.*w_bu,2) + w_n;

            
            % NOTE!
            % This is done with full matrices at the moment. 
            % Needs to be rewritten.
            Knn = B2'*B2;
            Knf = B2'*B;
            C = -L*L';
            for i=1:length(ind)
              La2 = gp_trcov(gp, xt(tstind{i},:), predcf) - B2(:,tstind{i})'*B2(:,tstind{i});
              Knn(ind{i},ind{i}) =  Knn(ind{i},ind{i}) + La2;
              Laa = gp_cov(gp, xt(tstind{i},:), x(ind{i},:),predcf) - B2(:,tstind{i})'*B(:,ind{i});
              Knf(tstind{i},ind{i}) =  Knf(tstind{i},ind{i}) + Laa;
              C(ind{i},ind{i}) = C(ind{i},ind{i}) + inv(La{i});
            end
            
            Covf = Knn - Knf * C * Knf';
            
            predcov = chol(Covf,'lower');
            Ef = repmat(Ef,1,nsamp);
            sampft = Ef + predcov*randn(size(Ef));
        end
        
        % ---------------------------
      case 'CS+FIC'
        n = size(x,1);
        n2 = size(xt,1);
        u = gp.X_u;
        m = length(u);

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
        
        K_fu = gp_cov(gp,x,u,predcf1);
        K_uu = gp_trcov(gp,u,predcf1);
        K_uu = (K_uu+K_uu')./2;
        K_nu=gp_cov(gp,xt,u,predcf1);

        Kcs_nf = gp_cov(gp, xt, x, predcf2);
        Kcs_nn = gp_trcov(gp, xt, predcf2);
        
        
        switch gp.latent_method
          case 'Laplace'
            
            
            [e, edata, eprior, f, L, a, La2] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
            

            deriv = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
            ntest=size(xt,1);
            
            % Calculate the predictive mean according to the type of
            % covariance functions used for making the prediction
            if ptype == 1
              Ef = K_nu*(K_uu\(K_fu'*deriv));
            elseif ptype == 2
              Ef = Kcs_nf*deriv;
            else 
              Ef = K_nu*(K_uu\(K_fu'*deriv)) + Kcs_nf*deriv;        
            end
            
            % evaluate also Lav if the prediction is made for training set
            if ~isempty(tstind)
              [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf1);
              Luu = chol(K_uu)';
              B=Luu\(K_fu');
              Qv_ff=sum(B.^2)';
              Lav = Kv_ff-Qv_ff;
            end
            
            % Add also Lav if the prediction is made for training set
            % and non-CS covariance function is used for prediction
            if ~isempty(tstind) && (ptype == 1 || ptype == 3)
              Ef(tstind) = Ef(tstind) + Lav.*deriv;
            end
            
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
            Covf = K - V'*V;
            
            predcov = chol(Covf,'lower');
            Ef = repmat(Ef,1,nsamp);
            sampft = Ef + predcov*randn(size(Ef));
            
          case 'EP'
            
            [e, edata, eprior, tautilde, nutilde, L, La, b] = gpep_e(gp_pak(gp), gp, x, y, 'z', z);

            p = b';
            ntest=size(xt,1);
            
            % Calculate the predictive mean according to the type of
            % covariance functions used for making the prediction
            if ptype == 1
              Ef = K_nu*(K_uu\(K_fu'*p));
            elseif ptype == 2
              Ef = Kcs_nf*p;
            else 
              Ef = K_nu*(K_uu\(K_fu'*p)) + Kcs_nf*p;        
            end

            % evaluate also Lav if the prediction is made for training set
            if ~isempty(tstind)
              [Kv_ff, Cv_ff] = gp_trvar(gp, xt(tstind,:), predcf1);
              Luu = chol(K_uu)';
              B=Luu\(K_fu');
              Qv_ff=sum(B.^2)';
              Lav = Kv_ff-Qv_ff;
            end
            
            % Add also Lav if the prediction is made for training set
            % and non-CS covariance function is used for prediction
            if ~isempty(tstind) && (ptype == 1 || ptype == 3)
              Ef(tstind) = Ef(tstind) + Lav.*p;
            end
            
            
            
            Luu = chol(K_uu)';
            B=Luu\(K_fu');   
            B2=Luu\(K_nu');   

            Knf = B2'*B + Kcs_nf;
            k = gp_trvar(gp,xt,predcf1);
            Knn = B2'*B2 + diag(k - sum(B2.^2)') + Kcs_nn;
            if ~isempty(tstind)
              Knf(tstind,:) = Knf(tstind,:) + diag(k(tstind) - sum(B.^2)');
            end
            
            Covf = Knn - Knf * ( inv(La) - L*L' ) * Knf';
            
            predcov = chol(Covf,'lower');
            Ef = repmat(Ef,1,nsamp);
            sampft = Ef + predcov*randn(size(Ef));
            
        end
    end
  end
elseif isstruct(gp) && numel(gp.jitterSigma2)>1
  % MCMC
  nmc=size(gp.jitterSigma2,1);
  % resample nsamp cases from nmc samples
  % deterministic resampling has low variance and small bias for
  % equal weights
  gi=resampdet(ones(nmc,1),nsamp,1);
  sampft=[];sampyt=[];
  for i1=1:nmc
    nsampi=sum(gi==i1);
    if nsampi>0
      Gp = take_nth(gp,i1);
      if isfield(Gp,'latent_method') && isequal(Gp.latent_method,'MCMC')
        Gp = rmfield(Gp,'latent_method');
      end
      if isfield(gp, 'latentValues') && ~isempty(gp.latentValues)
        % Non-Gaussian likelihood. The latent variables should be used in
        % place of observations
        y = gp.latentValues';
        ii=i1;
      else         
        ii=1;
      end
      if nargout<2
        tsampft = gp_rnd(Gp, x, y(:,ii), xt, 'nsamp', nsampi, ...
                         'z', z, 'zt', zt, 'predcf', predcf, ...
                         'tstind', tstind);
        sampft=[sampft tsampft];
      else
        [tsampft, tsampyt] = gp_rnd(Gp, x, y(:,ii), xt, 'nsamp', ...
                                    nsampi, 'z', z, 'zt', zt, ...
                                    'predcf', predcf, 'tstind', tstind);
        sampft=[sampft tsampft];
        sampyt=[sampyt tsampyt];
      end
    end
  end
elseif iscell(gp)
  % gp_ia
  ngp=length(gp);
  if isfield(gp{1},'ia_weight')
    for i1=1:length(gp)
      gw(i1)=gp{i1}.ia_weight;
    end
  else
    gw=ones(ngp,1);
  end
  % resample nsamp cases from nmc samples
  % strafied resampling has has higher variance than deterministic
  % resapmling, but has smaller bias, and thus it shoould be more
  % suitable for unequal weights
  gi=resampstr(gw,nsamp,1);
  sampft=[];sampyt=[];
  for i1 = 1:ngp
    nsampi=sum(gi==i1);
    if nsampi>0
      if nargout<2
        tsampft = gp_rnd(gp{i1}, x, y, xt, 'nsamp', nsampi, ...
                         'z', z, 'zt', zt, 'predcf', predcf, 'tstind', tstind);
        sampft=[sampft tsampft];
      else
        [tsampft, tsampyt] = gp_rnd(gp{i1}, x, y, xt, 'nsamp', nsampi, ...
                                    'z', z, 'zt', zt, 'predcf', predcf, 'tstind', tstind);
        sampft=[sampft tsampft];
        sampyt=[sampyt tsampyt];
      end
    end
  end
end
