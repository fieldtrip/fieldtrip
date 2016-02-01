function p_eff = gp_peff(gp, x, y, varargin);
%GP_PEFF  The effective number of parameters in GP model with focus 
%         on latent variables
%
%  Description
%    P_EFF = EP_PEFF(GP, X, Y) Takes the Gaussian process structure
%    GP, training inputs X and training outputs and returns the
%    effective number of parameters as defined by Spiegelhalter et
%    al. (2002).
%
%    NOTE! The effective number of parameters is evaluated with
%    focus on latent variable f. This means that the parameters th
%    (parameters of covariance function and likelihood) are
%    considered fixed. (See Spiegelhalter et al (2002) for
%    discussion on the parameters in focus in Bayesian model). 
%    Thus, the returned p_eff tells the effective number of latent
%    variables. This statistics is important for example when
%    assessing the goodness of Laplace or EP approximation in case
%    of non-Gaussian likelihood (See Vanhatalo et al for
%    discussion).
%
%    If you want to evaluate the effective number of parameters
%    with focus on parameters, see GP_DIC.
%
%    The effective number of parameters is approximated as follows:
%        p_eff = n - trace( K\C ),
%
%    where K is the prior covariance matrix and C the posterior
%    covariance matrix. This approximation is introduced by
%    Spiegelhalter et al. (2002) in equation (16). If the
%    likelihood is non-Gaussian and gp.latent_method is either
%    Laplace or EP, then C is the Laplace or EP approximation for
%    the posterior covariance.
%
%  References: 
%    Spiegelhalter, Best, Carlin and van der Linde (2002). 
%    Bayesian measures of model complexity and fit. J. R. 
%    Statist. Soc. B, 64(4):583-639.
%         
%    Vanhatalo, J., Pietil�inen V. and Vehtari, A. (2010). 
%    Approximate inference for disease mapping with sparse
%    Gaussian processes. Statistics in Medicine, 29(15):1580-1607.
%   
%  See also
%    GP_DIC, DEMO_MODELASSESMENT1
%   

% Copyright (c) 2009-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GP_PEFF';
  ip.addRequired('gp',@isstruct);
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.parse(gp, x, y, varargin{:});
  z=ip.Results.z;
  
  tn = size(x,1);

    
    if isfield(gp.lik.fh,'trcov')
      % a Gaussian likelihood
        
        switch gp.type
          case 'FULL'
                        
            [K, C] = gp_trcov(gp, x);
            L = chol(C);
            p_eff = trace( L\(L'\K) );
            
          case 'FIC'
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
            Luu = chol(K_uu)';
            
            % Evaluate Lambda (La) for specific model
            % Q_ff = K_fu*inv(K_uu)*K_fu'
            % Here we need only the diag(Q_ff), which is evaluated below
            B=Luu\(K_fu');
            Qv_ff=sum(B.^2)';
            Lav = Kv_ff-Qv_ff;
            Lav2 = Cv_ff-Qv_ff;

            iLaKfu = zeros(size(K_fu));  % f x u,
            n=size(x,1)
            for i=1:n
                iLaKfu(i,:) = K_fu(i,:)./Lav2(i);  % f x u
            end
            A = K_uu+K_fu'*iLaKfu;
            A = (A+A')./2;

            L = iLaKfu/chol(A);

            p_eff = sum(Lav./Lav2) + sum(sum( repmat(Lav2,1,m).\B'.*B' - L.*(L'*B'*B)' - L.*(L'.*repmat(Lav',m,1))', 2));
            
% $$$             % Check the result using full matrices
% $$$             C = B'*B + diag(Lav2);
% $$$             K = B'*B + diag(Lav);
% $$$             L = chol(C);
% $$$             p_eff = trace( L\(L'\K) );            
            
          case {'PIC' 'PIC_BLOCK'}
            u = gp.X_u;
            ind = gp.tr_index;
            if size(u,2) ~= size(x,2)
                u=u';
            end
            
            % Calculate some help matrices
            [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
            K_fu = gp_cov(gp, x, u);         % f x u
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
                La{i} = Kbl_ff - Qbl_ff;
                La2{i} = Cbl_ff - Qbl_ff;
                iLaKfu(ind{i},:) = La2{i}\K_fu(ind{i},:);    
            end
            A = K_uu+K_fu'*iLaKfu;
            A = (A+A')./2;            % Ensure symmetry
            L = iLaKfu/chol(A);
            
            p_eff = sum(sum(- L.*(L'*B'*B)',2));
            for i=1:length(ind)
                LLa2 = chol(La2{i});
                p_eff = p_eff + trace(LLa2\(LLa2'\La{i})) + trace( LLa2\(LLa2'\B(:,ind{i})'*B(:,ind{i})) - L(ind{i},:)*L(ind{i},:)'*La{i} );
            end
            
          case 'CS+FIC'
            u = gp.X_u;
            m = size(u,1);
            % Turn the inducing vector on right direction
            if size(u,2) ~= size(x,2)
                u=u';
            end

            % Indexes to all non-compact support and compact support covariances.
            cf1 = [];
            cf2 = [];
            
            ncf = length(gp.cf);
            % Loop through all covariance functions
            for i = 1:ncf        
                % Non-CS covariances
                if ~isfield(gp.cf{i},'cs') 
                    cf1 = [cf1 i];
                    % CS-covariances
                else
                    cf2 = [cf2 i];           
                end
            end

            [Kv_ff, Cv_ff] = gp_trvar(gp, x, cf1);  % f x 1  vector    
            K_fu = gp_cov(gp,x,u,cf1);         % f x u
            K_uu = gp_trcov(gp,u,cf1);    % u x u, noiseles covariance K_uu
            K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu
            
            Kcs = gp_trcov(gp, x, cf2);
            Luu = chol(K_uu)';
            B=Luu\(K_fu');       % u x f
            Qv_ff=sum(B.^2)';
            Lav = Cv_ff-Qv_ff;   % f x 1, Vector of diagonal elements
            Lav2 = Kv_ff-Qv_ff;   % f x 1, Vector of diagonal elements

            La = sparse(1:tn,1:tn,Lav,tn,tn) + Kcs;
            La2 = sparse(1:tn,1:tn,Lav2,tn,tn) + Kcs;
    
            iLaKfu = La\K_fu;
            A = K_uu+K_fu'*iLaKfu;
            A = (A+A')./2;     % Ensure symmetry
            L = iLaKfu/chol(A);
            
            LLa2 = L'*La2;
            LaB = La\B';
            LBB = L'*B'*B;
                        
            VD = ldlchol(La);
            spiLa = spinv(VD,1);
                
            p_eff = sum( sum(LaB.*B',2) + sum(spiLa.*La2,2) - sum(L.*LBB',2) - sum(L.*LLa2',2) );
            
% $$$             
% $$$             C = B'*B + Kcs + diag(Lav);
% $$$             K = B'*B + Kcs + diag(Lav2);
% $$$             
% $$$             L = chol(C);
% $$$             p_eff = trace( L\(L'\K) );
        end    
        
    else
      % ============================
      % A non Gaussian likelihood
      % ============================
        
        switch gp.type
          case 'FULL'
            switch gp.latent_method
              case 'EP'
                [e, edata, eprior, tautilde, nutilde, L] = gpep_e(gp_pak(gp), gp, x, y, 'z', z);
                
                % The prior variance
                K=gp_trcov(gp,x);
               
                if all(tautilde > 0) && ~isequal(gp.latent_opt.optim_method, 'robust-EP')
                    sqrttautilde = sqrt(tautilde);
                    Stildesqroot = sparse(1:tn, 1:tn, sqrttautilde, tn, tn);
                                        

                    if issparse(L)
                        p_eff = trace(Stildesqroot*ldlsolve(L, Stildesqroot*K));
                    else
                        p_eff = trace(Stildesqroot*(L'\(L\Stildesqroot)*K));
                    end
                else
                    C = L*L';
                    L = chol(K);
                    Varf = tn - trace(L\(L'\C));
                end
                
              case 'Laplace'
                [e, edata, eprior, f, L] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
                
                W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
                
                % Evaluate the prior variance
                K = gp_trcov(gp,x);
                
                if W >= 0
                    if issparse(K) && issparse(L)

                        sqrtW = sparse(1:tn, 1:tn, sqrt(W), tn, tn);
                        sqrtWK = sqrtW*K;
                        p_eff = trace(sqrtW*ldlsolve(L,sqrtWK));
                    else
                        W = diag(W);
                        p_eff = trace(sqrt(W)*(L'\(L\(sqrt(W)*K))));
                    end
                else
                    C = L*L';
                    L = chol(K);
                    p_eff = tn - trace(L\(L'\C));
                end
            end       
            
          case 'FIC'
            u = gp.X_u;
            m = size(u,1);
            % Calculate some help matrices
            [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
            K_fu = gp_cov(gp, x, u);   % f x u
            K_uu = gp_trcov(gp, u);     % u x u, noiseles covariance K_uu
            Luu = chol(K_uu)';
            B=Luu\(K_fu');
            Qv_ff=sum(B.^2)';
            Lav = Kv_ff-Qv_ff;
            
            switch gp.latent_method
              case 'EP'
                [e, edata, eprior, tautilde, nutilde, L, La, b] = gpep_e(gp_pak(gp), gp, x, y, 'z', z);

                k = gp_trvar(gp,x);
                
                p_eff = sum( sum( (repmat(1./La,1,m).*B').*B',2) ) - sum(sum(L.*((L'*B')*B)',2));
                Lav = k - sum(B.^2)';
                p_eff = p_eff + sum(Lav./La) - sum(sum(L.*L,2).*Lav);

% $$$                 C = diag(1./La) - L*L';
% $$$                 K = B'*B + diag(k - sum(B.^2)');
% $$$                 
% $$$                 p_eff = trace(C*K);
                
              case 'Laplace'
                [e, edata, eprior, f, L, a, La2] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
                
                W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
                La = W.*Lav;
                Lahat = 1 + La;
                sqrtW = sqrt(W);
                B = (repmat(sqrtW,1,m).*K_fu);
                
                % Components for (I + W^(1/2)*(Qff + La2)*W^(1/2))^(-1) = Lahat^(-1) - L2*L2'
                B2 = repmat(Lahat,1,m).\B;
                A2 = K_uu + B'*B2; A2=(A2+A2)/2;
                L2 = B2/chol(A2);
                
                BB=Luu\(K_fu');
                BB2=B/Luu';
                
                p_eff = sum(W./Lahat.*Lav);
                p_eff = p_eff + sum(sqrtW .* (sum((repmat(Lahat,1,m).\BB2).*BB',2)...
                                             - sum(L2.*(L2'.*repmat(sqrtW'.*La2',m,1))',2)...
                                             - sum(L2.*(L2'*BB2*BB)',2)) );
                
% $$$                 K = BB'*BB + diag(Lav);
% $$$                 sW = diag(sqrt(W));
% $$$                 B = eye(size(K)) + sW*K*sW;
% $$$                 L = chol(B)';
% $$$                                         
% $$$                 W = diag(W);
% $$$                 p_eff = trace(sqrt(W)*(L'\(L\(sqrt(W)*K))));
            end
                
          case {'PIC' 'PIC_BLOCK'}
            u = gp.X_u;
            K_fu = gp_cov(gp, x, u);
            K_uu = gp_trcov(gp, u);
            K_uu = (K_uu+K_uu')./2;
                        
            m = size(u,1);
            ind = gp.tr_index;
            Luu = chol(K_uu)';
            B=Luu\(K_fu');

            switch gp.latent_method
              case 'EP'
            
                [e, edata, eprior, tautilde, nutilde, L, La, b] = gpep_e(gp_pak(gp), gp, x, y, 'z', z);
        
                p_eff = - sum(sum(L.*((L'*B')*B)',2));

                for i=1:length(ind)
                    La1 = gp_trcov(gp, x(ind{i},:)) - B(:,ind{i})'*B(:,ind{i});
                    p_eff = p_eff + trace(La{i}\B(:,ind{i})'*B(:,ind{i}));
                    p_eff = p_eff + trace(La{i}\La1);
                    p_eff = p_eff - trace(L(ind{i},:)*L(ind{i},:)'*La1);
                end
                
              case 'Laplace'               
                [e, edata, eprior, f, L, a, La2] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
                
                
                iKuuKuf = K_uu\K_fu';

                % Evaluate the variance

                W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
                sqrtW = sqrt(W);
                
                % Components for (I + W^(1/2)*(Qff + La2)*W^(1/2))^(-1) = Lahat^(-1) - L2*L2'
                for i=1:length(ind)
                    La{i} = diag(sqrtW(ind{i}))*La2{i}*diag(sqrtW(ind{i}));
                    Lahat{i} = eye(size(La{i})) + La{i};
                    LLahat{i} = chol(Lahat{i});
                end
                sKfu = (repmat(sqrt(W),1,m).*K_fu);
                for i=1:length(ind)
                    iLasKfu(ind{i},:) = Lahat{i}\sKfu(ind{i},:);
                end
                A2 = K_uu + sKfu'*iLasKfu; A2=(A2+A2)/2;
                L2 = iLasKfu/chol(A2);
                
                
                p_eff = -sum(sqrtW.*sum(L2.*(L2'*(repmat(sqrtW,1,m).*B')*B)',2));
                for i=1:length(ind)
                    dsqrtW = diag(sqrtW(ind{i}));
                   p_eff = p_eff + trace( dsqrtW*(LLahat{i}\(LLahat{i}'\(dsqrtW*La2{i})))); 
                   p_eff = p_eff + trace( dsqrtW*(LLahat{i}\(LLahat{i}'\(dsqrtW*B(:,ind{i})'*B(:,ind{i}))))); 
                   p_eff = p_eff - trace(dsqrtW*L2(ind{i},:)*L2(ind{i},:)'*dsqrtW*La2{i});
                end
                
               
% $$$                 K = B'*B;
% $$$                 C = -L2*L2';
% $$$                 for i=1:length(ind)
% $$$                    K(ind{i},ind{i}) =  K(ind{i},ind{i}) + La2{i};
% $$$                    C(ind{i},ind{i}) =  C(ind{i},ind{i}) + inv(Lahat{i});
% $$$                 end
% $$$                                
% $$$                 p_eff = trace (diag(sqrtW)*C*diag(sqrtW)*K) ;
                
            end

            
          
          case 'CS+FIC'
            u = gp.X_u;
            [n,nin]=size(x);
            % Indexes to all non-compact support and compact support covariances.
            cf1 = [];
            cf2 = [];
            
            ncf = length(gp.cf);
            % Loop through all covariance functions
            for i = 1:ncf        
                % Non-CS covariances
                if ~isfield(gp.cf{i},'cs') 
                    cf1 = [cf1 i];
                    % CS-covariances
                else
                    cf2 = [cf2 i];           
                end
            end
            
            K_fu = gp_cov(gp,x,u,cf1);         % f x u
            K_uu = gp_trcov(gp,u,cf1);    % u x u, noiseles covariance K_uu
            K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu
                
            Kcs = gp_trcov(gp, x, cf2);
            Luu = chol(K_uu)';
            B=Luu\(K_fu');

            
            switch gp.latent_method
              case 'EP'

                [e, edata, eprior, tautilde, nutilde, L, La, b] = gpep_e(gp_pak(gp), gp, x, y, 'z', z);
            
                k = gp_trvar(gp,x,cf1);
                Lav = k - sum(B.^2)';
                La1 = Kcs + sparse(1:n, 1:n, Lav, n, n);
                
                issparse(La)
                VD = ldlchol(La);
                siLa = spinv(VD,1);
                p_eff = sum(sum(ldlsolve(VD,B').*B')) + sum(sum(siLa.*La1)) - sum(sum(L.*((L'*B')*B)',2));
                p_eff = p_eff - sum(sum(L.*(L'*La1)',2));

% $$$                 C = inv(La) - L*L';
% $$$                 K = B'*B + diag(k - sum(B.^2)') + Kcs;
% $$$                 
% $$$                 p_eff = trace(C*K);
                
              case 'Laplace'
                [e, edata, eprior, f, L, a, La2] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
                                
                W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
                sqrtW = sparse(1:tn,1:tn,sqrt(W),tn,tn);
                Lahat = sparse(1:tn,1:tn,1,tn,tn) + sqrtW*La2*sqrtW;
                B = sqrtW*K_fu;

                % Components for (I + W^(1/2)*(Qff + La2)*W^(1/2))^(-1) = Lahat^(-1) - L2*L2'
                B2 = Lahat\B;
                A2 = K_uu + B'*B2; A2=(A2+A2)/2;
                L2 = B2/chol(A2);

                BB = Luu\(K_fu');
                BB2 = B/Luu';

                VD = ldlchol(Lahat);
                spiLahat = spinv(VD,1);
                
                p_eff = sum(sum(sqrtW*spiLahat*sqrtW.*La2,2)) + sum(sqrtW * (sum(ldlsolve(VD, BB2).*BB',2) - sum(L2.*(L2'*sqrtW*La2)',2) - sum(L2.*(L2'*BB2*BB)',2)) );
                                
            end

            
            
            
        end
    end
    
    
end