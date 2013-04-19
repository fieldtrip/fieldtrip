function [Eft, Varft, lpyt, Eyt, Varyt] = gpla_loopred(gp, x, y, varargin)
%GPLA_LOOPRED  Leave-one-out predictions with Laplace approximation
%
%  Description
%    [EFT, VARFT, LPYT, EYT, VARYT] = GPLA_LOOPRED(GP, X, Y, OPTIONS)
%    takes a Gaussian process structure GP together with a matrix X
%    of training inputs and vector Y of training targets, and
%    evaluates the leave-one-out predictive distribution at inputs
%    X and returns means EFT and variances VARFT of latent
%    variables, the logarithm of the predictive densities PYT, and
%    the predictive means EYT and variances VARYT of observations
%    at input locations X.
%
%    OPTIONS is optional parameter-value pair
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case. 
%
%    Laplace leave-one-out is approximated in linear response style
%    by expressing the solutions for LOO problem in terms of
%    solution for the full problem. The computationally cheap
%    solution can be obtained by making the assumption that the
%    difference between these two solution is small such that their
%    difference may be treated as an Taylor expansion truncated at
%    first order (Winther et al 2012, in progress).
%
%  See also
%    GP_LOOPRED, GP_PRED
  
% Copyright (c) 2011-2012  Aki Vehtari, Ville Tolvanen

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPLA_LOOPRED';
  ip.addRequired('gp', @(x) isstruct(x));
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('method', 'lrs', @(x) ismember(x, {'lrs' 'cavity' 'inla'}))
  ip.parse(gp, x, y, varargin{:});
  z=ip.Results.z;
  method = ip.Results.method;
  [tn,nin] = size(x);
  
  switch method

    case 'lrs'
      % Manfred Opper and Ole Winther (2000). Gaussian Processes for
      % Classification: Mean-Field Algorithms. In Neural
      % Computation, 12(11):2655-2684.
      %
      % Ole Winther et al (2012). Work in progress.

      % latent posterior
      [f, sigm2ii] = gpla_pred(gp, x, y, 'z', z, 'tstind', []);
  
      deriv = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
      La = 1./-gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
      % really large values don't contribute, but make variance
      % computation unstable. 2e15 approx 1/(2*eps)
      La = min(La,2e15);
      
      switch gp.type
        case 'FULL' 
          % FULL GP (and compact support GP)
          K = gp_trcov(gp,x);
          Varft=1./diag(inv(K+diag(La)))-La;
          
        case 'FIC' 
          % FIC
          % Use inverse lemma for FIC low rank covariance matrix approximation
          % Code adapated from gp_pred
          
          u = gp.X_u;
          m = size(u,1);
          % Turn the inducing vector on right direction
          if size(u,2) ~= size(x,2)
            u=u';
          end
          [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
          K_fu = gp_cov(gp, x, u);   % f x u
          K_uu = gp_trcov(gp, u);     % u x u, noiseles covariance K_uu
          Luu = chol(K_uu,'lower');
          B=Luu\(K_fu');
          Qv_ff=sum(B.^2)';
          % Add also La to the vector of diagonal elements
          Lav = Cv_ff-Qv_ff + La;   % 1 x f, Vector of diagonal elements
          
          % iLaKfu = diag(inv(Lav))*K_fu = inv(La)*K_fu
          iLaKfu = zeros(size(K_fu));  % f x u,
          n=size(x,1);
          for i=1:n
            iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
          end
          A = K_uu+K_fu'*iLaKfu;
          A = (A+A')./2;
          L = iLaKfu/chol(A);
          
          %Varft=1./diag(inv(K+diag(La)))-La;
          Varft=1./(1./Lav - sum(L.^2,2))-La;
          
        case {'PIC' 'PIC_BLOCK'}
          % PIC
          % Use inverse lemma for PIC low rank covariance matrix approximation
          % Code adapated from gp_pred (here Lab is same La in gp_pred)
          
          u = gp.X_u;
          ind = gp.tr_index;
          if size(u,2) ~= size(x,2)
            % Turn the inducing vector on right direction
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
            % Add also La to the diagonal
            Lab{i} = Cbl_ff - Qbl_ff + diag(La(ind{i}));
            iLaKfu(ind{i},:) = Lab{i}\K_fu(ind{i},:);    
          end
          A = K_uu+K_fu'*iLaKfu;
          A = (A+A')./2;            % Ensure symmetry
          L = iLaKfu/chol(A);

          % From this on evaluate the prediction
          % See Snelson and Ghahramani (2007) for details
          n=size(y,1);
          iCv=zeros(n,1);
          for i=1:length(ind)
            iCv(ind{i},:) = diag(inv(Lab{i}));
          end

          %Varft=1./diag(inv(K+diag(La)))-La;
          Varft=1./(iCv - sum(L.^2,2))-La;

        case 'CS+FIC' 
          % CS+FIC
          % Use inverse lemma for CS+FIC
          % Code adapated from gp_pred (Here Las is same as La in gp_pred)
          
          u = gp.X_u;
          if size(u,2) ~= size(x,2)
            % Turn the inducing vector on right direction
            u=u';
          end
          
          n = size(x,1);
          m = size(u,1);
          ncf = length(gp.cf);
          
          % Indexes to all non-compact support and compact support covariances.
          cf1 = [];
          cf2 = [];

          % Loop through all covariance functions
          for i1 = 1:ncf        
            if ~isfield(gp.cf{i1},'cs') 
              % Non-CS covariances
              cf1 = [cf1 i1];
            else
              % CS-covariances
              cf2 = [cf2 i1];           
            end
          end
          
          % First evaluate needed covariance matrices
          % v defines that parameter is a vector
          [Kv_ff, Cv_ff] = gp_trvar(gp, x, cf1); % f x 1  vector    
          K_fu = gp_cov(gp, x, u, cf1);          % f x u
          K_uu = gp_trcov(gp, u, cf1);    % u x u, noiseles covariance K_uu
          K_uu = (K_uu+K_uu')./2;         % ensure the symmetry of K_uu

          Luu  = chol(K_uu)';
          % Evaluate the Lambda (La)
          % Q_ff = K_fu*inv(K_uu)*K_fu'
          B=Luu\(K_fu');       % u x f
          Qv_ff=sum(B.^2)';
          % Add also La to the vector of diagonal elements
          Lav = Cv_ff-Qv_ff + La;   % f x 1, Vector of diagonal elements

          K_cs = gp_trcov(gp,x,cf2);
          Las = sparse(1:n,1:n,Lav,n,n) + K_cs;

          iLaKfu = Las\K_fu;
          A = K_uu+K_fu'*iLaKfu;
          A = (A+A')./2;     % Ensure symmetry
          L = iLaKfu/chol(A);
          
          %Varft=1./diag(inv(K+diag(La)))-La;
          Varft=1./(diag(inv(Las)) - sum(L.^2,2))-La;
          
        otherwise
          error('Unknown type of Gaussian process')
      end
      
      % check if cavity variances are negative
      ii=find(Varft<0);
      if ~isempty(ii)
        warning('gpla_loopred: some LOO latent variances are negative');
        Varft(ii) = gp.jitterSigma2;
      end
      Eft=f-Varft.*deriv;

      if nargout==3
        lpyt = gp.lik.fh.predy(gp.lik, Eft, Varft, y, z);
      elseif nargout>3
        [lpyt,Eyt,Varyt] = gp.lik.fh.predy(gp.lik, Eft, Varft, y, z);
      end
      
    case 'cavity'
      % using EP equations

      % latent posterior
      [f, sigm2ii] = gpla_pred(gp, x, y, 'z', z, 'tstind', []);
      
      % "site parameters"
      W        = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
      deriv    = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
      sigm2_t  = 1./W;
      mu_t     = f + sigm2_t.*deriv;
      
      % "cavity parameters"
      sigma2_i = 1./(1./sigm2ii-1./sigm2_t);
      myy_i    = sigma2_i.*(f./sigm2ii-mu_t./sigm2_t);
      % check if cavity varianes are negative
      ii=find(sigma2_i<0);
      if ~isempty(ii)
        warning('gpla_loopred: some cavity variances are negative');
        sigma2_i(ii) = sigm2ii(ii);
        myy_i(ii) = f(ii);
      end
      
      % leave-one-out predictions
      Eft=myy_i;
      Varft=sigma2_i;

      if nargout==3
        lpyt = gp.lik.fh.predy(gp.lik, Eft, Varft, y, z);
      elseif nargout>3
        [lpyt,Eyt,Varyt] = gp.lik.fh.predy(gp.lik, Eft, Varft, y, z);
      end
      
    case 'inla'
      % Leonhard Held and Birgit Schr�dle and H�vard Rue (2010)
      % Posterior and Cross-validatory Predictive Checks: A
      % Comparison of MCMC and INLA. In (eds) Thomas Kneib and
      % Gerhard Tutz, Statistical Modelling and Regression
      % Structures, pp. 91-110. Springer.
      
      % latent posterior
      [f, sigm2ii, lp] = gpla_pred(gp, x, y, 'z', z, 'tstind', []);
      
      Eft = zeros(tn,1);
      Varft = zeros(tn,1);
      lpyt = zeros(tn,1);
      minf = f-6.*sqrt(sigm2ii);
      maxf = f+6.*sqrt(sigm2ii);
      for i=1:tn
        if isempty(z)
          z1 = [];
        else
          z1 = z(i);
        end
        [m0, m1, m2] = quad_moments(@(x) norm_pdf(x, f(i), sqrt(sigm2ii(i)))./llvec(gp.lik,y(i),x,z1), minf(i), maxf(i));
        Eft(i) = m1;
        Varft(i) = m2-Eft(i)^2;
        lpyt(i) = -log(m0);
      end

      if nargout>3
        [tmp,Eyt,Varyt] = gp.lik.fh.predy(gp.lik, Eft, Varft, y, z);
      end
      if sum((abs(lpyt)./abs(lp) > 5) == 1) > 0.1*tn;
        warning('Very bad predictive densities, gpla_loopred might not be reliable, check results!');
      end
      
  end

end

function expll = llvec(gplik, y, f, z)
  for i=1:size(f,2)
    expll(i) = exp(gplik.fh.ll(gplik, y, f(i), z));
  end
end
