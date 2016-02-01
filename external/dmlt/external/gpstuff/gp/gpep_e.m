function [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta] = gpep_e(w, gp, varargin)
%GPEP_E  Do Expectation propagation and return marginal log posterior estimate
%
%  Description
%    E = GPEP_E(W, GP, X, Y, OPTIONS) takes a GP structure GP
%    together with a matrix X of input vectors and a matrix Y of
%    target vectors, and finds the EP approximation for the
%    conditional posterior p(Y | X, th), where th is the
%    parameters. Returns the energy at th (see below). Each row of
%    X corresponds to one input vector and each row of Y
%    corresponds to one target vector.
%
%    [E, EDATA, EPRIOR] = GPEP_E(W, GP, X, Y, OPTIONS) returns also
%    the data and prior components of the total energy.
%
%    The energy is minus log posterior cost function for th:
%      E = EDATA + EPRIOR 
%        = - log p(Y|X, th) - log p(th),
%      where th represents the parameters (lengthScale,
%      magnSigma2...), X is inputs and Y is observations.
%
%    OPTIONS is optional parameter-value pair
%      z - optional observed quantity in triplet (x_i,y_i,z_i)
%          Some likelihoods may use this. For example, in case of
%          Poisson likelihood we have z_i=E_i, that is, expected
%          value for ith case.
%
%  See also
%    GP_SET, GP_E, GPEP_G, GPEP_PRED

%  Description 2
%    Additional properties meant only for internal use.
%  
%    GP = GPEP_E('init', GP) takes a GP structure GP and
%    initializes required fields for the EP algorithm.
%
%    GPEP_E('clearcache', GP) takes a GP structure GP and cleares
%    the internal cache stored in the nested function workspace
%
%    [e, edata, eprior, site_tau, site_nu, L, La2, b, muvec_i, sigm2vec_i]
%      = GPEP_E(w, gp, x, y, options)
%    returns many useful quantities produced by EP algorithm.
%
  
% Copyright (c) 2007  Jaakko Riihim�ki
% Copyright (c) 2007-2010  Jarno Vanhatalo
% Copyright (c) 2010 Heikki Peura
% Copyright (c) 2010-2012 Aki Vehtari
% Copyright (c) 2011 Pasi Jyl�nki

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
  
% parse inputs
  ip=inputParser;
  ip.FunctionName = 'GPEP_E';
  ip.addRequired('w', @(x) ...
                 isempty(x) || ...
                 (ischar(x) && ismember(x, {'init' 'clearcache'})) || ...
                 (isvector(x) && isreal(x) && all(isfinite(x))) || ...
                 all(isnan(x)));
  ip.addRequired('gp',@isstruct);
  ip.addOptional('x', [], @(x) isnumeric(x) && isreal(x) && all(isfinite(x(:))))
  ip.addOptional('y', [], @(x) isnumeric(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isnumeric(x) && isreal(x) && all(isfinite(x(:))))
  ip.parse(w, gp, varargin{:});
  x=ip.Results.x;
  y=ip.Results.y;
  z=ip.Results.z;
  
  if strcmp(w, 'init')
    % intialize cache
    ch = [];
    
    % return function handle to the nested function ep_algorithm
    % this way each gp has its own peristent memory for EP
    gp.fh.ne = @ep_algorithm;
    % set other function handles
    gp.fh.e=@gpep_e;
    gp.fh.g=@gpep_g;
    gp.fh.pred=@gpep_pred;
    gp.fh.jpred=@gpep_jpred;
    gp.fh.looe=@gpep_looe;
    gp.fh.loog=@gpep_loog;
    gp.fh.loopred=@gpep_loopred;
    e = gp;
    % remove clutter from the nested workspace
    clear w gp varargin ip x y z
  elseif strcmp(w, 'clearcache')
    % clear the cache
    gp.fh.ne('clearcache');
  else
    % call ep_algorithm using the function handle to the nested function
    % this way each gp has its own peristent memory for EP
    [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta] = gp.fh.ne(w, gp, x, y, z);
  end

  function [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta] = ep_algorithm(w, gp, x, y, z)
    
    if strcmp(w, 'clearcache')
      ch=[];
      return
    end
    
    switch gp.latent_opt.optim_method
      
      case 'basic-EP'
        
        % check whether saved values can be used
        if isempty(z)
          datahash=hash_sha512([x y]);
        else
          datahash=hash_sha512([x y z]);
        end
        if ~isempty(ch) && all(size(w)==size(ch.w)) && all(abs(w-ch.w)<1e-8) && isequal(datahash,ch.datahash)
          % The covariance function parameters or data haven't changed
          % so we can return the energy and the site parameters that are saved
          e = ch.e;
          edata = ch.edata;
          eprior = ch.eprior;
          tautilde = ch.tautilde;
          nutilde = ch.nutilde;
          L = ch.L;
          La2 = ch.La2;
          b = ch.b;
          muvec_i = ch.muvec_i;
          sigm2vec_i = ch.sigm2vec_i;
          logZ_i = ch.logZ_i;
          eta = ch.eta;
        else
          % The parameters or data have changed since
          % the last call for gpep_e. In this case we need to
          % re-evaluate the EP approximation
          gp=gp_unpak(gp, w);
          ncf = length(gp.cf);
          n = size(x,1);
          
          % EP iteration parameters
          iter=1;
          maxiter = gp.latent_opt.maxiter;
          tol = gp.latent_opt.tol;
          df = gp.latent_opt.df;
          nutilde = zeros(size(y));
          tautilde = zeros(size(y));
          muvec_i=zeros(size(y));
          sigm2vec_i=zeros(size(y));
          logZep_old=0; logZep=Inf;
          if ~isfield(gp,'meanf')
            mf = zeros(size(y));
          else
            [H,b_m,B_m]=mean_prep(gp,x,[]);
            mf = H'*b_m;
          end
          
          logM0 = zeros(n,1);
          muhat = zeros(n,1);
          sigm2hat = zeros(n,1);
          
          % =================================================
          % First Evaluate the data contribution to the error
          switch gp.type
            % ============================================================
            % FULL
            % ============================================================
            case 'FULL'   % A full GP
              [K,C] = gp_trcov(gp, x);
              if ~issparse(C)
                % The EP algorithm for full support covariance function
                if ~isfield(gp,'meanf')
                  Sigm = C;
                  meanfp=false;
                else
                  Sigm = C + H'*B_m*H;
                  meanfp=true;
                end
                
                % The EP -algorithm
                convergence=false;
                while iter<=maxiter && ~convergence
                  logZep_old=logZep;
                  logM0_old=logM0;
                  
                  if isequal(gp.latent_opt.init_prev, 'on') && iter==1 && ~isempty(ch) && all(size(w)==size(ch.w)) && all(abs(w-ch.w)<1) && isequal(datahash,ch.datahash)
                    tautilde=ch.tautilde;
                    nutilde=ch.nutilde;
                  else
                    if isequal(gp.latent_opt.parallel,'on')
                      % parallel-EP
                      % compute marginal and cavity parameters
                      dSigm=diag(Sigm);
                      tau=1./dSigm-tautilde;
                      nu = 1./dSigm.*mf-nutilde;
                      muvec_i=nu./tau;
                      sigm2vec_i=1./tau;
                      
                      % compute moments of tilted distributions
                      [logM0, muhat, sigm2hat] = gp.lik.fh.tiltedMoments(gp.lik, y, 1:n, sigm2vec_i, muvec_i, z);
                      if any(isnan(logM0))
                        [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                        return
                      end
                      
                      % update site parameters
                      deltatautilde=1./sigm2hat-tau-tautilde;
                      tautilde=tautilde+df.*deltatautilde;
                      deltanutilde=1./sigm2hat.*muhat-nu-nutilde;
                      nutilde=nutilde+df.*deltanutilde;
                    else
                      % sequential-EP
                      muvec_i = zeros(n,1); sigm2vec_i = zeros(n,1);
                      for i1=1:n
                        % Algorithm utilizing Cholesky updates
                        % This is numerically more stable but slower
                        % $$$                             % approximate cavity parameters
                        % $$$                             S11 = sum(Ls(:,i1).^2);
                        % $$$                             S1 = Ls'*Ls(:,i1);
                        % $$$                             tau_i=S11^-1-tautilde(i1);
                        % $$$                             nu_i=S11^-1*mf(i1)-nutilde(i1);
                        % $$$
                        % $$$                             mu_i=nu_i/tau_i;
                        % $$$                             sigm2_i=tau_i^-1;
                        % $$$
                        % $$$                             if sigm2_i < 0
                        % $$$                                 [ii i1]
                        % $$$                             end
                        % $$$
                        % $$$                             % marginal moments
                        % $$$                             [M0(i1), muhat, sigm2hat] = feval(gp.lik.fh.tiltedMoments, gp.lik, y, i1, sigm2_i, mu_i, z);
                        % $$$
                        % $$$                             % update site parameters
                        % $$$                             deltatautilde = sigm2hat^-1-tau_i-tautilde(i1);
                        % $$$                             tautilde(i1) = tautilde(i1)+deltatautilde;
                        % $$$                             nutilde(i1) = sigm2hat^-1*muhat-nu_i;
                        % $$$
                        % $$$                             upfact = 1./(deltatautilde^-1+S11);
                        % $$$                             if upfact > 0
                        % $$$                                 Ls = cholupdate(Ls, S1.*sqrt(upfact), '-');
                        % $$$                             else
                        % $$$                                 Ls = cholupdate(Ls, S1.*sqrt(-upfact));
                        % $$$                             end
                        % $$$                             Sigm = Ls'*Ls;
                        % $$$                             mf=Sigm*nutilde;
                        % $$$
                        % $$$                             muvec_i(i1,1)=mu_i;
                        % $$$                             sigm2vec_i(i1,1)=sigm2_i;
                        
                        % Algorithm as in Rasmussen and Williams 2006
                        % approximate cavity parameters
                        Sigmi=Sigm(:,i1);
                        Sigmii=Sigmi(i1);
                        tau_i=1/Sigmii-tautilde(i1);
                        nu_i = 1/Sigmii*mf(i1)-nutilde(i1);
                        mu_i=nu_i/tau_i;
                        sigm2_i=1/tau_i;
                        
                        % marginal moments
                        [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2_i, mu_i, z);
                        if isnan(logM0(i1))
                          [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                          return
                        end
                        % update site parameters
                        deltatautilde=sigm2hat(i1)^-1-tau_i-tautilde(i1);
                        tautilde(i1)=tautilde(i1)+df*deltatautilde;
                        deltanutilde=sigm2hat(i1)^-1*muhat(i1)-nu_i-nutilde(i1);
                        nutilde(i1)=nutilde(i1)+df*deltanutilde;
                        
                        % Update mean and variance after each site update (standard EP)
                        ds = deltatautilde/(1+deltatautilde*Sigmii);
                        Sigm = Sigm - ((ds*Sigmi)*Sigmi');
                        %Sigm = Sigm - ((ds*Sigm(:,i1))*Sigm(:,i1)');
                        % The below is how Rasmussen and Williams
                        % (2006) do the update. The above version is
                        % more robust.
                        %ds = deltatautilde^-1+Sigm(i1,i1);
                        %ds = (Sigm(:,i1)/ds)*Sigm(:,i1)';
                        %Sigm = Sigm - ds;
                        %Sigm=Sigm-(deltatautilde^-1+Sigm(i1,i1))^-1*(Sigm(:,i1)*Sigm(:,i1)');
                        
                        if ~meanfp
                          mf=Sigm*nutilde;
                        else
                          mf=Sigm*(C\(H'*b_m)+nutilde);
                        end
                        
                        muvec_i(i1)=mu_i;
                        sigm2vec_i(i1)=sigm2_i;
                      end
                    end
                  end
                  
                  % Recompute the approximate posterior parameters
                  % parallel- and sequential-EP
                  Stilde=tautilde;
                  Stildesqr=sqrt(Stilde);
                  
                  if ~meanfp % zero mean function used
                    % NOTICE! upper triangle matrix! cf. to
                    % line 13 in the algorithm 3.5, p. 58.
                    
                    %B=eye(n)+Stildesqr*C*Stildesqr;
                    B=bsxfun(@times,bsxfun(@times,Stildesqr,C),Stildesqr');
                    B(1:n+1:end)=B(1:n+1:end)+1;
                    [L,notpositivedefinite] = chol(B,'lower');
                    if notpositivedefinite
                      [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                      return
                    end
                    %V=(L\Stildesqr)*C;
                    V=L\bsxfun(@times,Stildesqr,C);
                    Sigm=C-V'*V; 
                    mf=Sigm*nutilde;

                    % Compute the marginal likelihood
                    % Direct formula (3.65):
                    % Sigmtilde=diag(1./tautilde);
                    % mutilde=inv(Stilde)*nutilde;
                    %
                    % logZep=-0.5*log(det(Sigmtilde+K))-0.5*mutilde'*inv(K+Sigmtilde)*mutilde+
                    %         sum(log(normcdf(y.*muvec_i./sqrt(1+sigm2vec_i))))+
                    %         0.5*sum(log(sigm2vec_i+1./tautilde))+
                    %         sum((muvec_i-mutilde).^2./(2*(sigm2vec_i+1./tautilde)))
                    
                    % 4. term & 1. term
                    term41=0.5*sum(log(1+tautilde.*sigm2vec_i))-sum(log(diag(L)));
                    
                    % 5. term (1/2 element) & 2. term
                    T=1./sigm2vec_i;
                    Cnutilde = C*nutilde;
                    L2 = V*nutilde;
                    term52 = nutilde'*Cnutilde - L2'*L2 - (nutilde'./(T+Stilde)')*nutilde;
                    term52 = term52.*0.5;
                    
                    % 5. term (2/2 element)
                    term5=0.5*muvec_i'.*(T./(Stilde+T))'*(Stilde.*muvec_i-2*nutilde);
                    
                    % 3. term
                    term3 = sum(logM0);
                    
                    logZep = -(term41+term52+term5+term3);
                    iter=iter+1;
                    
                  else
                    % mean function used
                    % help variables
                    hBh = H'*B_m*H;
                    C_t = C + hBh;
                    CHb  = C\H'*b_m;
                    S   = diag(Stildesqr.^2);
                    %B = eye(n)+Stildesqroot*C*Stildesqroot;
                    B=bsxfun(@times,bsxfun(@times,Stildesqr,C),Stildesqr');
                    B(1:n+1:end)=B(1:n+1:end)+1;
                    %B_h = eye(n) + Stildesqroot*C_t*Stildesqroot;
                    B_h=bsxfun(@times,bsxfun(@times,Stildesqr,C_t),Stildesqr');
                    B_h(1:n+1:end)=B_h(1:n+1:end)+1;
                    % L to return, without the hBh term
                    [L,notpositivedefinite]=chol(B,'lower');
                    if notpositivedefinite
                      [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                      return
                    end
                    % L for the calculation with mean term
                    [L_m,notpositivedefinite]=chol(B_h,'lower');
                    if notpositivedefinite
                      [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                      return
                    end
                    
                    % Recompute the approximate posterior parameters
                    % parallel- and sequential-EP
                    
                    %V=(L_m\Stildesqroot)*C_t;
                    V=L_m\bsxfun(@times,Stildesqr,C_t);
                    Sigm=C_t-V'*V;
                    mf=Sigm*(CHb+nutilde);
                    
                    T=1./sigm2vec_i;
                    Cnutilde = (C_t - S^-1)*(S*H'*b_m-nutilde);
                    L2 = V*(S*H'*b_m-nutilde);
                    
                    Stildesqroot = diag(Stildesqr);
                    zz   = Stildesqroot*(L'\(L\(Stildesqroot*C)));
                    % inv(K + S^-1)*S^-1
                    Ks  = eye(size(zz)) - zz;
                    
                    % 5. term (1/2 element)
                    term5_1  = 0.5.*((nutilde'*S^-1)./(T.^-1+Stilde.^-1)')*(S^-1*nutilde);
                    % 2. term
                    term2    = 0.5.*((S*H'*b_m-nutilde)'*Cnutilde - L2'*L2);
                    % 4. term
                    term4    = 0.5*sum(log(1+tautilde.*sigm2vec_i));
                    % 1. term
                    term1    = -1.*sum(log(diag(L_m)));
                    % 3. term
                    term3    = sum(logM0);
                    % 5. term (2/2 element)
                    term5    = 0.5*muvec_i'.*(T./(Stilde+T))'*(Stilde.*muvec_i-2*nutilde);

                    logZep = -(term4+term1+term5_1+term5+term2+term3);
                    
                    iter=iter+1;
                    
                  end
                  convergence=max(abs(logM0_old-logM0))<tol && abs(logZep_old-logZep)<tol;
                end

              else
                % EP algorithm for compactly supported covariance function
                % (C is a sparse matrix)
                p = analyze(K);
                r(p) = 1:n;
                if ~isempty(z)
                  z = z(p,:);
                end
                y = y(p);
                K = K(p,p);
                
                Inn = sparse(1:n,1:n,1,n,n);
                sqrtS = sparse(1:n,1:n,0,n,n);
                mf = zeros(size(y));
                sigm2 = zeros(size(y));
                dSigm=full(diag(K));
                gamma = zeros(size(y));
                VD = sparse(1:n,1:n,1,n,n);
                
                % The EP -algorithm
                convergence=false;
                while iter<=maxiter && ~convergence
                  logZep_old=logZep;
                  logM0_old=logM0;
                  
                  if isequal(gp.latent_opt.parallel,'on')
                    % parallel-EP
                    % approximate cavity parameters
                    sqrtSK = ssmult(sqrtS, K);
                    tttt = ldlsolve(VD,sqrtSK);
                    sigm2 = full(diag(K) - sum(sqrtSK.*tttt)');
                    mf = gamma - tttt'*sqrtS*gamma;
                    tau=1./sigm2-tautilde;
                    nu = 1./sigm2.*mf-nutilde;
                    muvec_i=nu./tau;
                    sigm2vec_i=1./tau;
                    % compute moments of tilted distributions
                    for i1=1:n
                      [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2vec_i(i1), muvec_i(i1), z);
                    end
                    if any(isnan(logM0))
                      [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                      return
                    end
                    % update site parameters
                    deltatautilde=1./sigm2hat-tau-tautilde;
                    tautilde=tautilde+df*deltatautilde;
                    deltanutilde=1./sigm2hat.*muhat-nu-nutilde;
                    nutilde=nutilde+df*deltanutilde;
                    gamma = gamma + sum(bsxfun(@times,K,df.*deltanutilde'),2);
                  else
                    % sequential-EP
                    muvec_i = zeros(n,1); sigm2vec_i = zeros(n,1);
                    for i1=1:n
                      % approximate cavity parameters
                      Ki1 = K(:,i1);
                      sqrtSKi1 = ssmult(sqrtS, Ki1);
                      tttt = ldlsolve(VD,sqrtSKi1);
                      sigm2(i1) = Ki1(i1) - sqrtSKi1'*tttt;
                      mf(i1) = gamma(i1) - tttt'*sqrtS*gamma;
                      
                      tau_i=sigm2(i1)^-1-tautilde(i1);
                      nu_i=sigm2(i1)^-1*mf(i1)-nutilde(i1);
                      
                      mu_i=nu_i/tau_i;
                      sigm2_i=tau_i^-1;
                      
                      % marginal moments
                      [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2_i, mu_i, z);
                      
                      % update site parameters
                      tautilde_old = tautilde(i1);
                      deltatautilde=sigm2hat(i1)^-1-tau_i-tautilde(i1);
                      tautilde(i1)=tautilde(i1)+df*deltatautilde;
                      deltanutilde=sigm2hat(i1)^-1*muhat(i1)-nu_i-nutilde(i1);
                      nutilde(i1)=nutilde(i1)+df*deltanutilde;
                      gamma = gamma + Ki1.*df*deltanutilde;
                      
                      % Update the LDL decomposition
                      sqrtS(i1,i1) = sqrt(tautilde(i1));
                      sqrtSKi1(i1) = sqrt(tautilde(i1)).*Ki1(i1);
                      D2_n = sqrtSKi1.*sqrtS(i1,i1) + Inn(:,i1);
                      
                      if tautilde_old == 0
                        VD = ldlrowupdate(i1,VD,VD(:,i1),'-');
                        VD = ldlrowupdate(i1,VD,D2_n,'+');
                      else
                        VD = ldlrowmodify(VD, D2_n, i1);
                      end
                      
                      muvec_i(i1,1)=mu_i;
                      sigm2vec_i(i1,1)=sigm2_i;
                    end

                  end
                  
                  % Recompute the approximate posterior parameters
                  % parallel- and sequential-EP
                  sqrtS = sparse(1:n,1:n,sqrt(tautilde),n,n);
                  KsqrtS = ssmult(K,sqrtS);
                  B = ssmult(sqrtS,KsqrtS) + Inn;
                  [VD, notpositivedefinite] = ldlchol(B);
                  if notpositivedefinite
                    [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                    return
                  end
                  Knutilde = K*nutilde;
                  mf = Knutilde - KsqrtS*ldlsolve(VD,sqrtS*Knutilde);
                  
                  % Compute the marginal likelihood
                  % 4. term & 1. term
                  term41=0.5*sum(log(1+tautilde.*sigm2vec_i)) - 0.5.*sum(log(diag(VD)));
                  
                  % 5. term (1/2 element) & 2. term
                  T=1./sigm2vec_i;
                  term52 = nutilde'*mf - (nutilde'./(T+tautilde)')*nutilde;
                  term52 = term52.*0.5;
                  
                  % 5. term (2/2 element)
                  term5=0.5*muvec_i'.*(T./(tautilde+T))'*(tautilde.*muvec_i-2*nutilde);
                  
                  % 3. term
                  term3 = sum(logM0);
                  
                  logZep = -(term41+term52+term5+term3);
                  
                  iter=iter+1;
                  
                convergence=max(abs(logM0_old-logM0))<tol && abs(logZep_old-logZep)<tol;
                %[iter-1 max(abs(muhat-mf)./abs(mf)) max(abs(sqrt(sigm2hat)-s)./abs(s)) max(abs(logM0_old-logM0)) abs(logZep_old-logZep)]
                %[iter-1 max(abs(muhat-mf)./abs(mf)) max(abs(logM0_old-logM0)) abs(logZep_old-logZep)]
                end
                % Reorder all the returned and stored values
                
                B = B(r,r);
                nutilde = nutilde(r);
                tautilde = tautilde(r);
                muvec_i = muvec_i(r);
                sigm2vec_i = sigm2vec_i(r);
                logM0 = logM0(r);
                mf = mf(r);
                y = y(r);
                if ~isempty(z)
                  z = z(r,:);
                end
                [L, notpositivedefinite] = ldlchol(B);
                if notpositivedefinite
                  [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                  return
                end
                
              end
              edata = logZep;
              % Set something into La2
              La2 = B;
              b = 0;
              
              % ============================================================
              % FIC
              % ============================================================
            case 'FIC'
              u = gp.X_u;
              m = size(u,1);
              
              % First evaluate needed covariance matrices
              % v defines that parameter is a vector
              [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % f x 1  vector
              K_fu = gp_cov(gp, x, u);           % f x u
              K_uu = gp_trcov(gp, u);     % u x u, noiseles covariance K_uu
              K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu
              [Luu, notpositivedefinite] = chol(K_uu, 'lower');
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              % Evaluate the Lambda (La)
              % Q_ff = K_fu*inv(K_uu)*K_fu'
              % Here we need only the diag(Q_ff), which is evaluated below
              B=Luu\(K_fu');       % u x f
              Qv_ff=sum(B.^2)';
              Lav = Cv_ff-Qv_ff;   % f x 1, Vector of diagonal elements
              % iLaKfu = diag(iLav)*K_fu = inv(La)*K_fu
              % First some helper parameters
              iLaKfu = zeros(size(K_fu));  % f x u,
              for i=1:n
                iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
              end
              A = K_uu+K_fu'*iLaKfu;  A = (A+A')./2;     % Ensure symmetry
              [A, notpositivedefinite] = chol(A);
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              L = iLaKfu/A;
              Lahat = 1./Lav;
              I = eye(size(K_uu));
              
              [R0, notpositivedefinite] = chol(inv(K_uu));
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              R = R0;
              P = K_fu;
              mf = zeros(size(y));
              eta = zeros(size(y));
              gamma = zeros(size(K_uu,1),1);
              D_vec = Lav;
              Ann=0;
              
              % The EP -algorithm
              convergence=false;
              while iter<=maxiter && ~convergence
                logZep_old=logZep;
                logM0_old=logM0;
                  
                if isequal(gp.latent_opt.parallel,'on')
                  % parallel-EP
                  % approximate cavity parameters
                  Ann = D_vec+sum((P*R').^2,2);
                  mf = eta + sum(bsxfun(@times,P,gamma'),2);
                  tau = 1./Ann-tautilde;
                  nu = 1./Ann.*mf-nutilde;
                  muvec_i=nu./tau;
                  sigm2vec_i=1./tau;
                  % compute moments of tilted distributions
                  for i1=1:n
                    [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2vec_i(i1), muvec_i(i1), z);
                  end
                  if any(isnan(logM0))
                    [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                    return
                  end
                  % update site parameters
                  deltatautilde=1./sigm2hat-tau-tautilde;
                  tautilde=tautilde+df*deltatautilde;
                  deltanutilde=1./sigm2hat.*muhat-nu-nutilde;
                  nutilde=nutilde+df*deltanutilde;
                else
                  % sequential-EP
                  muvec_i = zeros(n,1); sigm2vec_i = zeros(n,1);
                  for i1=1:n
                    % approximate cavity parameters
                    pn = P(i1,:)';
                    Ann = D_vec(i1) + sum((R*pn).^2);
                    tau_i = Ann^-1-tautilde(i1);
                    mf(i1) = eta(i1) + pn'*gamma;
                    nu_i = Ann^-1*mf(i1)-nutilde(i1);
                    
                    mu_i=nu_i/tau_i;
                    sigm2_i=tau_i^-1;
                  
                    % marginal moments
                    [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2_i, mu_i, z);
                  
                    % update site parameters
                    deltatautilde = sigm2hat(i1)^-1-tau_i-tautilde(i1);
                    tautilde(i1) = tautilde(i1)+df*deltatautilde;
                    deltanutilde = sigm2hat(i1)^-1*muhat(i1)-nu_i - nutilde(i1);
                    nutilde(i1) = nutilde(i1)+df*deltanutilde;
                  
                    % Update the parameters
                    dn = D_vec(i1);
                    D_vec(i1) = D_vec(i1) - deltatautilde.*D_vec(i1).^2 ./ (1+deltatautilde.*D_vec(i1));
                    P(i1,:) = pn' - (deltatautilde.*dn ./ (1+deltatautilde.*dn)).*pn';
                    updfact = deltatautilde./(1 + deltatautilde.*Ann);
                    if updfact > 0
                      RtRpnU = R'*(R*pn).*sqrt(updfact);
                      R = cholupdate(R, RtRpnU, '-');
                    elseif updfact < 0
                      RtRpnU = R'*(R*pn).*sqrt(abs(updfact));
                      R = cholupdate(R, RtRpnU, '+');
                    end
                    eta(i1) = eta(i1) + (deltanutilde - deltatautilde.*eta(i1)).*dn./(1+deltatautilde.*dn);
                    gamma = gamma + (deltanutilde - deltatautilde.*mf(i1))./(1+deltatautilde.*dn) * R'*(R*pn);
                    %                            mf = eta + P*gamma;
                  
                    % Store cavity parameters
                    muvec_i(i1,1)=mu_i;
                    sigm2vec_i(i1,1)=sigm2_i;
                  end
                end
                
                % Recompute the approximate posterior parameters
                % parallel- and sequential-EP
                temp1 = (1+Lav.*tautilde).^(-1);
                D_vec = temp1.*Lav;
                R0P0t = R0*K_fu';
                temp2 = zeros(size(R0P0t));
%                for i2 = 1:length(temp1)
%                  P(i2,:) = temp1(i2).*K_fu(i2,:);
%                  temp2(:,i2) = R0P0t(:,i2).*tautilde(i2).*temp1(i2);
%                end
%                R = chol(inv(eye(size(R0)) + temp2*R0P0t')) * R0;
                P=bsxfun(@times,temp1,K_fu);
                temp2=bsxfun(@times,(tautilde.*temp1)',R0P0t);
                temp2=temp2*R0P0t';
                temp2(1:m+1:end)=temp2(1:m+1:end)+1;
                R = chol(inv(temp2)) * R0;
                eta = D_vec.*nutilde;
                gamma = R'*(R*(P'*nutilde));
                mf = eta + P*gamma;
                
                % Compute the marginal likelihood, see FULL model for
                % details about equations
                Lahat = 1./Lav + tautilde;
                Lhat = bsxfun(@rdivide,L,Lahat);
                H = I-L'*Lhat;
                B = H\L';
                Bhat = B./repmat(Lahat',m,1);
                
                % 4. term & 1. term
                Stildesqroot=sqrt(tautilde);
                D = Stildesqroot.*Lav.*Stildesqroot + 1;
                SsqrtKfu = K_fu.*repmat(Stildesqroot,1,m);
                AA = K_uu + (SsqrtKfu'./repmat(D',m,1))*SsqrtKfu; AA = (AA+AA')/2;
                [AA, notpositivedefinite] = chol(AA,'lower');
                if notpositivedefinite
                  [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                  return
                end
                term41 = - 0.5*sum(log(1+tautilde.*sigm2vec_i)) - sum(log(diag(Luu))) + sum(log(diag(AA))) + 0.5.*sum(log(D));
                
                % 5. term (1/2 element) & 2. term
                T=1./sigm2vec_i;
                term52 = -0.5*( (nutilde./Lahat)'*nutilde + (nutilde'*Lhat)*(Bhat*nutilde) - (nutilde./(T+tautilde))'*nutilde);
                
                % 5. term (2/2 element)
                term5 = - 0.5*muvec_i'.*(T./(tautilde+T))'*(tautilde.*muvec_i-2*nutilde);
                
                % 3. term
                term3 = -sum(logM0);
                
                logZep = term41+term52+term5+term3;
                
                iter=iter+1;
                convergence=max(abs(logM0_old-logM0))<tol && abs(logZep_old-logZep)<tol;
              end
              edata = logZep;
              %L = iLaKfu;
              
              % b'  = (La + Kfu*iKuu*Kuf + 1./S)*1./S * nutilde
              %     = (S - S * (iLa - L*L' + S)^(-1) * S) * 1./S
              %     = I - S * (Lahat - L*L')^(-1)
              % L   = S*Kfu * (Lav + 1./S)^(-1) / chol(K_uu + SsqrtKfu'*(Lav + 1./S)^(-1)*SsqrtKfu)
              % La2 = D./S = Lav + 1./S,
              %
              % The way evaluations are done is numerically more stable
              % See equations (3.71) and (3.72) in Rasmussen and Williams (2006)
              b = nutilde'.*(1 - Stildesqroot./Lahat.*Stildesqroot)' - (nutilde'*Lhat)*Bhat.*tautilde';    % part of eq. (3.71)
              L = ((repmat(Stildesqroot,1,m).*SsqrtKfu)./repmat(D',m,1)')/AA';                             % part of eq. (3.72)
              La2 = 1./(Stildesqroot./D.*Stildesqroot);                                                    % part of eq. (3.72)
              D = D_vec;
              
              % ============================================================
              % PIC
              % ============================================================
            case {'PIC' 'PIC_BLOCK'}
              ind = gp.tr_index;
              u = gp.X_u;
              m = length(u);
              
              % First evaluate needed covariance matrices
              % v defines that parameter is a vector
              K_fu = gp_cov(gp, x, u);         % f x u
              
              K_uu = gp_trcov(gp, u);    % u x u, noiseles covariance K_uu
              K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu
              [Luu, notpositivedefinite] = chol(K_uu, 'lower');
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              % Evaluate the Lambda (La)
              % Q_ff = K_fu*inv(K_uu)*K_fu'
              % Here we need only the diag(Q_ff), which is evaluated below
              B=Luu\(K_fu');       % u x f
              
              % First some helper parameters
              iLaKfu = zeros(size(K_fu));  % f x u
              for i=1:length(ind)
                Qbl_ff = B(:,ind{i})'*B(:,ind{i});
                [Kbl_ff, Cbl_ff] = gp_trcov(gp, x(ind{i},:));
                Labl{i} = Cbl_ff - Qbl_ff;
                [Llabl, notpositivedefinite] = chol(Labl{i});
                if notpositivedefinite
                  [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                  return
                end
                iLaKfu(ind{i},:) = Llabl\(Llabl'\K_fu(ind{i},:));
              end
              A = K_uu+K_fu'*iLaKfu;
              A = (A+A')./2;     % Ensure symmetry
              [A, notpositivedefinite] = chol(A);
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              L = iLaKfu/A;
              I = eye(size(K_uu));
              
              [R0, notpositivedefinite] = chol(inv(K_uu));
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              R = R0;
              P = K_fu;
              R0P0t = R0*K_fu';
              mf = zeros(size(y));
              eta = zeros(size(y));
              gamma = zeros(size(K_uu,1),1);
              D = Labl;
              Ann=0;
              
              % The EP -algorithm
              convergence=false;
              while iter<=maxiter && ~convergence
                logZep_old=logZep;
                logM0_old=logM0;
                
                if isequal(gp.latent_opt.parallel,'on')
                  % parallel-EP
                  % approximate cavity parameters
                  for bl=1:length(ind)
                    bl_ind = ind{bl};
                    Pbl=P(bl_ind,:);
                    Ann = diag(D{bl}) +sum((Pbl*R').^2,2);
                    tau(bl_ind,1) = 1./Ann-tautilde(bl_ind);
                    mf(bl_ind,1) = eta(bl_ind) + sum(bsxfun(@times,Pbl,gamma'),2);
                    nu(bl_ind,1) = 1./Ann.*mf(bl_ind)-nutilde(bl_ind);
                  end
                  muvec_i=nu./tau;
                  sigm2vec_i=1./tau;
                  % compute moments of tilted distributions
                  for i1=1:n
                    [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2vec_i(i1), muvec_i(i1), z);
                  end
                  if any(isnan(logM0))
                    [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                    return
                  end
                  % update site parameters
                  deltatautilde = 1./sigm2hat-tau-tautilde;
                  tautilde = tautilde+df*deltatautilde;
                  deltanutilde = 1./sigm2hat.*muhat-nu-nutilde;
                  nutilde = nutilde+df*deltanutilde;;
                else
                
                  muvec_i = zeros(n,1); sigm2vec_i = zeros(n,1);
                  for bl=1:length(ind)
                    bl_ind = ind{bl};
                    for in=1:length(bl_ind)
                      i1 = bl_ind(in);
                      % approximate cavity parameters
                      Dbl = D{bl}; dn = Dbl(in,in); pn = P(i1,:)';
                      Ann = dn + sum((R*pn).^2);
                      tau_i = Ann^-1-tautilde(i1);
                      mf(i1) = eta(i1) + pn'*gamma;
                      nu_i = Ann^-1*mf(i1)-nutilde(i1);
                      
                      mu_i=nu_i/tau_i;
                      sigm2_i=tau_i^-1;
                      
                      % marginal moments
                      [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2_i, mu_i, z);
                      
                      % update site parameters
                      deltatautilde = sigm2hat(i1)^-1-tau_i-tautilde(i1);
                      tautilde(i1) = tautilde(i1)+df*deltatautilde;
                      deltanutilde = sigm2hat(i1)^-1*muhat(i1)-nu_i - nutilde(i1);
                      nutilde(i1) = nutilde(i1) + df*deltanutilde;
                    
                      % Update the parameters
                      Dblin = Dbl(:,in);
                      Dbl = Dbl - deltatautilde ./ (1+deltatautilde.*dn) * Dblin*Dblin';
                      %Dbl = inv(inv(Dbl) + diag(tautilde(bl_ind)));
                      P(bl_ind,:) = P(bl_ind,:) - ((deltatautilde ./ (1+deltatautilde.*dn)).* Dblin)*pn';
                      updfact = deltatautilde./(1 + deltatautilde.*Ann);
                      if updfact > 0
                        RtRpnU = R'*(R*pn).*sqrt(updfact);
                        R = cholupdate(R, RtRpnU, '-');
                      elseif updfact < 0
                        RtRpnU = R'*(R*pn).*sqrt(abs(updfact));
                        R = cholupdate(R, RtRpnU, '+');
                      end
                      eta(bl_ind) = eta(bl_ind) + (deltanutilde - deltatautilde.*eta(i1))./(1+deltatautilde.*dn).*Dblin;
                      gamma = gamma + (deltanutilde - deltatautilde.*mf(i1))./(1+deltatautilde.*dn) * (R'*(R*pn));
                      %mf = eta + P*gamma;
                        
                      D{bl} = Dbl;
                      % Store cavity parameters
                      muvec_i(i1,1)=mu_i;
                      sigm2vec_i(i1,1)=sigm2_i;
                    end
                  end
                end
                
                % Recompute the approximate posterior parameters
                % parallel- and sequential-EP
                temp2 = zeros(size(R0P0t));
                
                Stildesqroot=sqrt(tautilde);
                for i=1:length(ind)
                  sdtautilde = diag(Stildesqroot(ind{i}));
                  Dhat = sdtautilde*Labl{i}*sdtautilde + eye(size(Labl{i}));
                  [Ldhat{i}, notpositivedefinite] = chol(Dhat);
                  if notpositivedefinite
                    [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                    return
                  end
                  D{i} = Labl{i} - Labl{i}*sdtautilde*(Ldhat{i}\(Ldhat{i}'\sdtautilde*Labl{i}));
                  P(ind{i},:) = D{i}*(Labl{i}\K_fu(ind{i},:));
                  
                  temp2(:,ind{i}) = R0P0t(:,ind{i})*sdtautilde/Dhat*sdtautilde;
                  eta(ind{i}) = D{i}*nutilde(ind{i});
                end
                R = chol(inv(eye(size(R0)) + temp2*R0P0t')) * R0;
                gamma = R'*(R*(P'*nutilde));
                mf = eta + P*gamma;
                
                % Compute the marginal likelihood, see FULL model for
                % details about equations
                %
                % First some helper parameters
                for i = 1:length(ind)
                  Lhat(ind{i},:) = D{i}*L(ind{i},:);
                end
                H = I-L'*Lhat;
                B = H\L';
                
                % Compute the marginal likelihood, see FULL model for
                % details about equations
                term41 = 0; term52 = 0;
                for i=1:length(ind)
                  Bhat(:,ind{i}) = B(:,ind{i})*D{i};
                  SsqrtKfu(ind{i},:) = bsxfun(@times,K_fu(ind{i},:),Stildesqroot(ind{i}));
                  %SsqrtKfu(ind{i},:) = gtimes(K_fu(ind{i},:),Stildesqroot(ind{i}));
                  iDSsqrtKfu(ind{i},:) = Ldhat{i}\(Ldhat{i}'\SsqrtKfu(ind{i},:));
                  term41 = term41 + sum(log(diag(Ldhat{i})));
                  term52 = term52 + nutilde(ind{i})'*(D{i}*nutilde(ind{i}));
                  
                end
                AA = K_uu + SsqrtKfu'*iDSsqrtKfu; AA = (AA+AA')/2;
                [AA, notpositivedefinite] = chol(AA,'lower');
                if notpositivedefinite
                  [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                  return
                end
                term41 = term41 - 0.5*sum(log(1+tautilde.*sigm2vec_i)) - sum(log(diag(Luu))) + sum(log(diag(AA)));
                
                % 5. term (1/2 element) & 2. term
                T=1./sigm2vec_i;
                term52 = -0.5*( term52 + (nutilde'*Lhat)*(Bhat*nutilde) - (nutilde./(T+tautilde))'*nutilde);
                
                % 5. term (2/2 element)
                term5 = - 0.5*muvec_i'.*(T./(tautilde+T))'*(tautilde.*muvec_i-2*nutilde);
                
                % 3. term
                term3 = -sum(logM0);
                
                logZep = term41+term52+term5+term3;
                
                iter=iter+1;
                convergence=max(abs(logM0_old-logM0))<tol && abs(logZep_old-logZep)<tol;
              end
              edata = logZep;
              
              b = zeros(1,n);
              for i=1:length(ind)
                b(ind{i}) = nutilde(ind{i})'*D{i};
                La2{i} = inv(diag(Stildesqroot(ind{i}))*(Ldhat{i}\(Ldhat{i}'\diag(Stildesqroot(ind{i})))));
              end
              b = nutilde' - ((b + (nutilde'*Lhat)*Bhat).*tautilde');
              
              L = (repmat(Stildesqroot,1,m).*iDSsqrtKfu)/AA';
              
              % ============================================================
              % CS+FIC
              % ============================================================
            case 'CS+FIC'
              u = gp.X_u;
              m = length(u);
              cf_orig = gp.cf;
              
              cf1 = {};
              cf2 = {};
              j = 1;
              k = 1;
              for i = 1:ncf
                if ~isfield(gp.cf{i},'cs')
                  cf1{j} = gp.cf{i};
                  j = j + 1;
                else
                  cf2{k} = gp.cf{i};
                  k = k + 1;
                end
              end
              gp.cf = cf1;
              
              % First evaluate needed covariance matrices
              % v defines that parameter is a vector
              [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % f x 1  vector
              K_fu = gp_cov(gp, x, u);           % f x u
              K_uu = gp_trcov(gp, u);            % u x u, noiseles covariance K_uu
              K_uu = (K_uu+K_uu')./2;            % ensure the symmetry of K_uu
              [Luu, notpositivedefinite] = chol(K_uu, 'lower');
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              
              % Evaluate the Lambda (La)
              % Q_ff = K_fu*inv(K_uu)*K_fu'
              B=Luu\(K_fu');       % u x f
              Qv_ff=sum(B.^2)';
              Lav = Cv_ff-Qv_ff;   % f x 1, Vector of diagonal elements
              
              gp.cf = cf2;
              K_cs = gp_trcov(gp,x);
              La = sparse(1:n,1:n,Lav,n,n) + K_cs;
              gp.cf = cf_orig;
              
              % clear unnecessary variables
              clear K_cs; clear Qv_ff; clear Kv_ff; clear Cv_ff; clear Lav;
              
              % Find fill reducing permutation and permute all the
              % matrices
              p = analyze(La);
              r(p) = 1:n;
              if ~isempty(z)
                z = z(p,:);
              end
              y = y(p);
              La = La(p,p);
              K_fu = K_fu(p,:);
              
              [VD, notpositivedefinite] = ldlchol(La);
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              iLaKfu = ldlsolve(VD,K_fu);
              A = K_uu+K_fu'*iLaKfu; A = (A+A')./2;     % Ensure symmetry
              [A, notpositivedefinite] = chol(A);
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              L = iLaKfu/A;
              
              I = eye(size(K_uu));
              
              Inn = sparse(1:n,1:n,1,n,n);
              sqrtS = sparse(1:n,1:n,0,n,n);
              [R0, notpositivedefinite] = chol(inv(K_uu));
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              R = R0;
              P = K_fu;
              R0P0t = R0*K_fu';
              mf = zeros(size(y));
              eta = zeros(size(y));
              gamma = zeros(size(K_uu,1),1);
              Ann=0;
              LasqrtS = La*sqrtS;
              [VD, notpositivedefinite] = ldlchol(Inn);
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              
              % The EP -algorithm
              convergence=false;
              while iter<=maxiter && ~convergence
                logZep_old=logZep;
                logM0_old=logM0;
                
                if isequal(gp.latent_opt.parallel,'on')
                  % parallel-EP
                  % approximate cavity parameters
                  tttt = ldlsolve(VD,ssmult(sqrtS,La));
                  D_vec = full(diag(La) - sum(LasqrtS'.*tttt)');
                  Ann = D_vec+sum((P*R').^2,2);
                  mf = eta + sum(bsxfun(@times,P,gamma'),2);
                  tau = 1./Ann-tautilde;
                  nu = 1./Ann.*mf-nutilde;
                  muvec_i=nu./tau;
                  sigm2vec_i= 1./tau;
                  % compute moments of tilted distributions
                  for i1=1:n
                    [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2vec_i(i1), muvec_i(i1), z);
                  end
                  if any(isnan(logM0))
                    [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                    return
                  end
                  % update site parameters
                  deltatautilde=1./sigm2hat-tau-tautilde;
                  tautilde=tautilde+df*deltatautilde;
                  deltanutilde=1./sigm2hat.*muhat-nu-nutilde;
                  nutilde=nutilde+df*deltanutilde;
                else
                  % sequential-EP
                  muvec_i = zeros(n,1); sigm2vec_i = zeros(n,1);
                  for i1=1:n
                    % approximate cavity parameters
                    tttt = ldlsolve(VD,ssmult(sqrtS,La(:,i1)));
                    Di1 =  La(:,i1) - ssmult(LasqrtS,tttt);
                    
                    dn = Di1(i1);
                    pn = P(i1,:)';
                    Ann = dn + sum((R*pn).^2);
                    tau_i = Ann^-1-tautilde(i1);
                    mf(i1) = eta(i1) + pn'*gamma;
                    nu_i = Ann^-1*mf(i1)-nutilde(i1);
                    
                    mu_i=nu_i/tau_i;
                    sigm2_i= tau_i^-1;  % 1./tau_i;  %
                    
                    % marginal moments
                    [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2_i, mu_i, z);
                    
                    % update site parameters
                    deltatautilde = sigm2hat(i1)^-1-tau_i-tautilde(i1);
                    tautilde(i1) = tautilde(i1)+df*deltatautilde;
                    deltanutilde = sigm2hat(i1)^-1*muhat(i1)-nu_i - nutilde(i1);
                    nutilde(i1) = nutilde(i1) + df*deltanutilde;
                  
                    % Update the parameters
                    P = P - ((deltatautilde ./ (1+deltatautilde.*dn)).* Di1)*pn';
                    updfact = deltatautilde./(1 + deltatautilde.*Ann);
                    if updfact > 0
                      RtRpnU = R'*(R*pn).*sqrt(updfact);
                      R = cholupdate(R, RtRpnU, '-');
                    elseif updfact < 0
                      RtRpnU = R'*(R*pn).*sqrt(abs(updfact));
                      R = cholupdate(R, RtRpnU, '+');
                    end
                    eta = eta + (deltanutilde - deltatautilde.*eta(i1))./(1+deltatautilde.*dn).*Di1;
                    gamma = gamma + (deltanutilde - deltatautilde.*mf(i1))./(1+deltatautilde.*dn) * (R'*(R*pn));
                    
                    % Store cavity parameters
                    muvec_i(i1,1)=mu_i;
                    sigm2vec_i(i1,1)=sigm2_i;
                    
                    D2_o = ssmult(sqrtS,LasqrtS(:,i1)) + Inn(:,i1);
                    sqrtS(i1,i1) = sqrt(tautilde(i1));
                    LasqrtS(:,i1) = La(:,i1).*sqrtS(i1,i1);
                    D2_n = ssmult(sqrtS,LasqrtS(:,i1)) + Inn(:,i1);
                    
                    if tautilde(i1) - deltatautilde == 0
                      VD = ldlrowupdate(i1,VD,VD(:,i1),'-');
                      VD = ldlrowupdate(i1,VD,D2_n,'+');
                    else
                      VD = ldlrowmodify(VD, D2_n, i1);
                    end
                  end
                end
                
                % Recompute the approximate posterior parameters
                % parallel- and sequential-EP
                sqrtS = sparse(1:n,1:n,sqrt(tautilde),n,n);
                sqrtSLa = ssmult(sqrtS,La);
                D2 = ssmult(sqrtSLa,sqrtS) + Inn;
                LasqrtS = ssmult(La,sqrtS);
                [VD, notpositivedefinite] = ldlchol(D2);
                if notpositivedefinite
                  [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                  return
                end
                
                SsqrtKfu = sqrtS*K_fu;
                iDSsqrtKfu = ldlsolve(VD,SsqrtKfu);
                P = K_fu - sqrtSLa'*iDSsqrtKfu;
                R = chol(inv( eye(size(R0)) + R0P0t*sqrtS*ldlsolve(VD,sqrtS*R0P0t'))) * R0;
                eta = La*nutilde - sqrtSLa'*ldlsolve(VD,sqrtSLa*nutilde);
                gamma = R'*(R*(P'*nutilde));
                mf = eta + P*gamma;
                
                % Compute the marginal likelihood,
                Lhat = La*L - sqrtSLa'*ldlsolve(VD,sqrtSLa*L);
                H = I-L'*Lhat;
                B = H\L';
                Bhat = B*La - ldlsolve(VD,sqrtSLa*B')'*sqrtSLa;
                
                % 4. term & 1. term
                AA = K_uu + SsqrtKfu'*iDSsqrtKfu; AA = (AA+AA')/2;
                [AA, notpositivedefinite] = chol(AA,'lower');
                if notpositivedefinite
                  [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                  return
                end
                term41 = - 0.5*sum(log(1+tautilde.*sigm2vec_i)) - sum(log(diag(Luu))) + sum(log(diag(AA))) + 0.5*sum(log(diag(VD)));
                
                % 5. term (1/2 element) & 2. term
                T=1./sigm2vec_i;
                term52 = -0.5*( nutilde'*(eta) + (nutilde'*Lhat)*(Bhat*nutilde) - (nutilde./(T+tautilde))'*nutilde);
                
                % 5. term (2/2 element)
                term5 = - 0.5*muvec_i'.*(T./(tautilde+T))'*(tautilde.*muvec_i-2*nutilde);
                
                % 3. term
                term3 = -sum(logM0);
                
                logZep = term41+term52+term5+term3;
                
                iter=iter+1;
                convergence=max(abs(logM0_old-logM0))<tol && abs(logZep_old-logZep)<tol;
              end
              edata = logZep;
              
              % b'  = (K_fu/K_uu*K_fu' + La + diag(1./tautilde)) \ (tautilde.\nutilde)
              % L   = S*Kfu * (Lav + 1./S)^(-1) / chol(K_uu + SsqrtKfu'*(Lav + 1./S)^(-1)*SsqrtKfu)
              % La2 = D./S = Lav + 1./S,
              %
              % The way evaluations are done is numerically more stable than with inversion of S (tautilde)
              % See equations (3.71) and (3.72) in Rasmussen and Williams (2006)
              b = nutilde' - ((eta' + (nutilde'*Lhat)*Bhat).*tautilde');
              
              L = (sqrtS*iDSsqrtKfu)/AA';
              La2 = sqrtS\D2/sqrtS;
              
              % Reorder all the returned and stored values
              b = b(r);
              L = L(r,:);
              La2 = La2(r,r);
              D = La(r,r);
              nutilde = nutilde(r);
              tautilde = tautilde(r);
              logM0 = logM0(r);
              muvec_i = muvec_i(r);
              sigm2vec_i = sigm2vec_i(r);
              mf = mf(r);
              P = P(r,:);
              y = y(r);
              if ~isempty(z)
                z = z(r,:);
              end
              % ============================================================
              % DTC,VAR
              % ============================================================
            case {'DTC' 'VAR' 'SOR'}
              % First evaluate needed covariance matrices
              % v defines that parameter is a vector
              u = gp.X_u;
              m = size(u,1);
              
              % First evaluate needed covariance matrices
              % v defines that parameter is a vector
              [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % f x 1  vector
              K_fu = gp_cov(gp, x, u);           % f x u
              K_uu = gp_trcov(gp, u);     % u x u, noiseles covariance K_uu
              K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu
              [Luu, notpositivedefinite] = chol(K_uu, 'lower');
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              % Evaluate the Lambda (La)
              % Q_ff = K_fu*inv(K_uu)*K_fu'
              % Here we need only the diag(Q_ff), which is evaluated below
              B=Luu\(K_fu');       % u x f
              
              
              Phi = B';
              m = size(Phi,2);
              
              R = eye(m,m);
              P = Phi;
              mf = zeros(size(y));
              gamma = zeros(m,1);
              Ann=0;
              
              % The EP -algorithm
              convergence=false;
              while iter<=maxiter && ~convergence
                logZep_old=logZep;
                logM0_old=logM0;
                
                if isequal(gp.latent_opt.parallel,'on')
                  % parallel-EP
                  % approximate cavity parameters
                  Ann = sum((P*R').^2,2);
                  mf = sum(bsxfun(@times,Phi,gamma'),2);%phi'*gamma;
                  tau = 1./Ann-tautilde;
                  nu = 1./Ann.*mf-nutilde;
                  muvec_i=nu./tau;
                  sigm2vec_i= 1./tau;
                  % compute moments of tilted distributions
                  for i1=1:n
                    [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2vec_i(i1), muvec_i(i1), z);
                  end
                  if any(isnan(logM0))
                    [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                    return
                  end
                  % update site parameters
                  deltatautilde=1./sigm2hat-tau-tautilde;
                  tautilde=tautilde+df*deltatautilde;
                  deltanutilde=1./sigm2hat.*muhat-nu-nutilde;
                  nutilde=nutilde+df*deltanutilde;
                else
                  % sequential-EP
                  muvec_i = zeros(n,1); sigm2vec_i = zeros(n,1);
                  for i1=1:n
                    % approximate cavity parameters
                    phi = Phi(i1,:)';
                    Ann = sum((R*phi).^2);
                    tau_i = Ann^-1-tautilde(i1);
                    mf(i1) = phi'*gamma;
                    nu_i = Ann^-1*mf(i1)-nutilde(i1);
                    
                    mu_i=nu_i/tau_i;
                    sigm2_i=tau_i^-1;
                    
                    % marginal moments
                    [logM0(i1), muhat(i1), sigm2hat(i1)] = gp.lik.fh.tiltedMoments(gp.lik, y, i1, sigm2_i, mu_i, z);
                    
                    % update site parameters
                    deltatautilde = sigm2hat(i1)^-1-tau_i-tautilde(i1);
                    tautilde(i1) = tautilde(i1)+df*deltatautilde;
                    deltanutilde = sigm2hat(i1)^-1*muhat(i1)-nu_i - nutilde(i1);
                    nutilde(i1) = nutilde(i1) + df*deltanutilde;
                  
                    % Update the parameters
                    lnn = sum((R*phi).^2);
                    updfact = deltatautilde/(1 + deltatautilde*lnn);
                    if updfact > 0
                      RtLphiU = R'*(R*phi).*sqrt(updfact);
                      R = cholupdate(R, RtLphiU, '-');
                    elseif updfact < 0
                      RtLphiU = R'*(R*phi).*sqrt(updfact);
                      R = cholupdate(R, RtLphiU, '+');
                    end
                    gamma = gamma - R'*(R*phi)*(deltatautilde*mf(i1)-deltanutilde);
                    % Store cavity parameters
                    muvec_i(i1,1)=mu_i;
                    sigm2vec_i(i1,1)=sigm2_i;
                  end
                end
                
                % Recompute the approximate posterior parameters
                % parallel- and sequential-EP
                R = chol(inv(eye(m,m) + Phi'*(repmat(tautilde,1,m).*Phi)));
                gamma = R'*(R*(Phi'*nutilde));
                mf = Phi*gamma;
                
                % Compute the marginal likelihood, see FULL model for
                % details about equations
                % 4. term & 1. term
                Stildesqroot=sqrt(tautilde);
                SsqrtPhi = Phi.*repmat(Stildesqroot,1,m);
                AA = eye(m,m) + SsqrtPhi'*SsqrtPhi; AA = (AA+AA')/2;
                [AA, notpositivedefinite] = chol(AA,'lower');
                if notpositivedefinite
                  [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                  return
                end
                term41 = - 0.5*sum(log(1+tautilde.*sigm2vec_i)) + sum(log(diag(AA)));
                
                % 5. term (1/2 element) & 2. term
                T=1./sigm2vec_i;
                bb = nutilde'*Phi;
                bb2 = bb*SsqrtPhi';
                bb3 = bb2*SsqrtPhi/AA';
                term52 = -0.5*( bb*bb' - bb2*bb2' + bb3*bb3' - (nutilde./(T+tautilde))'*nutilde);
                
                % 5. term (2/2 element)
                term5 = - 0.5*muvec_i'.*(T./(tautilde+T))'*(tautilde.*muvec_i-2*nutilde);
                
                % 3. term
                term3 = -sum(logM0);
                
                logZep = term41+term52+term5+term3;
                
                iter=iter+1;
                convergence=max(abs(logM0_old-logM0))<tol && abs(logZep_old-logZep)<tol;
              end
              edata = logZep;
              %L = iLaKfu;
              if strcmp(gp.type,'VAR')
                Qv_ff = sum(B.^2)';
                edata = edata + 0.5*sum((Kv_ff-Qv_ff).*tautilde);
              end
              
              temp = Phi*(SsqrtPhi'*(SsqrtPhi*bb'));
              %                b = Phi*bb' - temp + Phi*(SsqrtPhi'*(SsqrtPhi*(AA'\(AA\temp))));
              
              b = nutilde - bb2'.*Stildesqroot + repmat(tautilde,1,m).*Phi*(AA'\bb3');
              b = b';
              
              %                 StildeKfu = zeros(size(K_fu));  % f x u,
              %                 for i=1:n
              %                     StildeKfu(i,:) = K_fu(i,:).*tautilde(i);  % f x u
              %                 end
              %                 A = K_uu+K_fu'*StildeKfu;  A = (A+A')./2;     % Ensure symmetry
              %                 A = chol(A);
              %                 L = StildeKfu/A;
              L = repmat(tautilde,1,m).*Phi/AA';
              %L = repmat(tautilde,1,m).*K_fu/AA';
              mu=nutilde./tautilde;
              %b = nutilde - mu'*L*L'*mu;
              %b=b';
              La2 = 1./tautilde;
              D = 0;
              
            otherwise
              error('Unknown type of Gaussian process!')
          end
          
          % ==================================================
          % Evaluate the prior contribution to the error from
          % covariance functions and likelihood
          % ==================================================
          
          % Evaluate the prior contribution to the error from covariance
          % functions
          eprior = 0;
          for i=1:ncf
            gpcf = gp.cf{i};
            eprior = eprior - gpcf.fh.lp(gpcf);
          end
          
          % Evaluate the prior contribution to the error from likelihood
          % functions
          if isfield(gp.lik, 'p')
            lik = gp.lik;
            eprior = eprior - lik.fh.lp(lik);
          end
          
          % The last things to do
          if isfield(gp.latent_opt, 'display') && ismember(gp.latent_opt.display,{'final','iter'})
            fprintf('GPEP_E: Number of iterations in EP: %d \n', iter-1)
          end
          
          e = edata + eprior;
          logZ_i = logM0(:);
          eta = [];
          
          % store values to the cache
          ch.w = w;
          ch.e = e;
          ch.edata = edata;
          ch.eprior = eprior;
          ch.tautilde = tautilde;
          ch.nutilde = nutilde;
          ch.L = L;
          ch.La2 = La2;
          ch.b = b;
          ch.muvec_i = muvec_i;
          ch.sigm2vec_i = sigm2vec_i;
          ch.logZ_i = logZ_i;
          ch.eta = eta;
          ch.datahash=datahash;
          
          global iter_lkm
          iter_lkm=iter;
        end
        
      case 'robust-EP'
        
        %   function [e,edata,eprior,tau_q,nu_q,L, La2, b, muvec_i,sigm2vec_i,Z_i, eta] = ep_algorithm2(w,gp,x,y,z)
        %
        %     if strcmp(w, 'clearcache')
        %       ch=[];
        %       return
        %     end
        
        % check whether saved values can be used
        if isempty(z)
          datahash=hash_sha512([x y]);
        else
          datahash=hash_sha512([x y z]);
        end
        
        if ~isempty(ch) && all(size(w)==size(ch.w)) && all(abs(w-ch.w) < 1e-8) && isequal(datahash,ch.datahash)
          % The covariance function parameters haven't changed so just
          % return the Energy and the site parameters that are saved
          
          e = ch.e;
          edata = ch.edata;
          eprior = ch.eprior;
          L = ch.L;
          La2 = ch.La2;
          b = ch.b;
          nutilde = ch.nu_q;
          tautilde = ch.tau_q;
          eta = ch.eta;
          muvec_i = ch.muvec_i;
          sigm2vec_i = ch.sigm2vec_i;
          logZ_i = ch.logZ_i;
          
        else
          % The parameters or data have changed since
          % the last call for gpep_e. In this case we need to
          % re-evaluate the EP approximation
          
          % preparations
          ninit=gp.latent_opt.ninit; % max number of initial parallel iterations
          maxiter=gp.latent_opt.maxiter; % max number of double-loop iterations
          max_ninner=gp.latent_opt.max_ninner; % max number of inner loop iterations in the double-loop algorithm
          tolStop=gp.latent_opt.tolStop; % converge tolerance
          tolUpdate=gp.latent_opt.tolUpdate; % tolerance for the EP site updates
          tolInner=gp.latent_opt.tolInner; % inner loop energy tolerance
          tolGrad=gp.latent_opt.tolGrad; % minimum gradient (g) decrease in the search direction, abs(g_new)<tolGrad*abs(g)
          Vc_lim=gp.latent_opt.cavity_var_lim; % limit for the cavity variance Vc, Vc < Vc_lim*diag(K)
          df0=gp.latent_opt.df; % the intial damping factor
          eta1=gp.latent_opt.eta; % the initial fraction parameter
          eta2=gp.latent_opt.eta2; % the secondary fraction parameter          
          display=gp.latent_opt.display; % control the display
          
          gp=gp_unpak(gp,w);
          likelih=gp.lik;
          ncf = length(gp.cf);
          n=length(y);      
          pvis=0;
          
          eta=repmat(eta1,n,1);  % the initial vector of fraction parameters
          fh_tm=@(si,m_c,V_c,eta) likelih.fh.tiltedMoments2(likelih,y,si,V_c,m_c,z,eta);
          
          switch gp.type
            case 'FULL'
              % prior covariance
              K = gp_trcov(gp, x);
              
            case 'FIC'
              % Sparse
              u = gp.X_u;
              m = size(u,1);
              K_uu = gp_trcov(gp,u);
              K_uu = (K_uu + K_uu')./2;
              K_fu = gp_cov(gp,x,u);
              [Kv_ff, Cv_ff] = gp_trvar(gp,x);
              [Luu, notpositivedefinite] = chol(K_uu, 'lower');
              if notpositivedefinite
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
              B=Luu\(K_fu');
              Qv_ff=sum(B.^2)';
              Sf = [];
              Sf2 = [];
              L2 = [];
          end
          
          % prior (zero) initialization
          [nu_q,tau_q]=deal(zeros(n,1));
          
          % initialize the q-distribution (the multivariate Gaussian posterior approximation)
          switch gp.type
            case 'FULL'
              [mf,Sf,lnZ_q]=evaluate_q(nu_q,tau_q,K,display);
              Vf = diag(Sf);
            case 'FIC'
              [mf,Vf,lnZ_q]=evaluate_q2(nu_q,tau_q,Luu, K_fu, Kv_ff, Qv_ff, display);
            otherwise
              error('Robust-EP not implemented for this type of GP!');
          end
          
          % initialize the surrogate distribution (the independent Gaussian marginal approximations)
          nu_s=mf./Vf;
          tau_s=1./Vf;
          lnZ_s=0.5*sum( (-log(tau_s) +nu_s.^2 ./tau_s)./eta ); % minus 0.5*log(2*pi)./eta
          
          % initialize r-distribution (the tilted distributions)
          [lnZ_r,lnZ_i,m_r,V_r]=evaluate_r(nu_q,tau_q,eta,fh_tm,nu_s,tau_s,display);
          
          % initial energy (lnZ_ep)
          e = lnZ_q + lnZ_r -lnZ_s;
          
          if ismember(display,{'final','iter'})
            fprintf('\nInitial energy: e=%.4f, hyperparameters:\n',e)
            fprintf('Cov:%s \n',sprintf(' %.2g,',gp_pak(gp,'covariance')))
            fprintf('Lik:%s \n',sprintf(' %.2g,',gp_pak(gp,'likelihood')))
          end
          
          if isfinite(e) % do not run the algorithm if the prior energy is not defined
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % initialize with ninit rounds of parallel EP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % EP search direction
            up_mode='ep'; % choose the moment matching
            [dnu_q,dtau_q]=ep_update_dir(mf,Vf,m_r,V_r,eta,up_mode,tolUpdate);
            
            convergence=false; % convergence indicator
            df=df0; % initial damping factor
            tol_m=zeros(1,2); % absolute moment tolerances
            switch gp.type
              case 'FULL'
                tauc_min=1./(Vc_lim*diag(K)); % minimum cavity precision
              case 'FIC'
                tauc_min=1./(Vc_lim*Cv_ff);
            end
            % Adjust damping by setting an upper limit (Vf_mult) to the increase
            % of the marginal variance
            Vf_mult=2;
            i1=0;
            while i1<ninit
              i1=i1+1;
              
              %%%%%%%%%%%%%%%%%%%
              % the damped update
              dfi=df(ones(n,1));
              temp=(1/Vf_mult-1)./Vf;
              ii2=df*dtau_q<temp;
              if any(ii2)
                dfi(ii2)=temp(ii2)./dtau_q(ii2);
              end
              
              % proposal site parameters
              nu_q2=nu_q+dfi.*dnu_q;
              tau_q2=tau_q+dfi.*dtau_q;
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%
              % a proposal q-distribution
              switch gp.type
                case 'FULL'
                  [mf2,Sf2,lnZ_q2,L1,L2]=evaluate_q(nu_q2,tau_q2,K,display);
                  Vf2 = diag(Sf2);
                case 'FIC'
                  [mf2,Vf2,lnZ_q2,L1,L2]=evaluate_q2(nu_q2,tau_q2,Luu, K_fu, Kv_ff, Qv_ff, display);
                otherwise
                  error('Robust-EP not implemented for this type of GP!');
              end
              
              % check that the new cavity variances do not exceed the limit
              tau_s2=1./Vf2;
              pcavity=all( (tau_s2-eta.*tau_q2 )>=tauc_min);
              if isempty(L2) || ~pcavity
                % In case of too small cavity precisions, half the step size
                df=df*0.5;
                if df<0.1,
                  % If mediocre damping is not sufficient, proceed to
                  % the double-loop algorithm
                  break
                else
                  if ismember(display,{'iter'})
                    fprintf('%d, e=%.6f, dm=%.4f, dV=%.4f, increasing damping to df=%g.\n',i1,e,tol_m(1),tol_m(2),df)
                  end
                  continue
                end
              end
              
              % a proposal surrogate distribution
              nu_s2=mf2./Vf2;
              lnZ_s2=0.5*sum( (-log(tau_s2) +nu_s2.^2 ./tau_s2)./eta );
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%
              % a proposal r-distribution
              [lnZ_r2,lnZ_i2,m_r2,V_r2,p]=evaluate_r(nu_q2,tau_q2,eta,fh_tm,nu_s2,tau_s2,display);
              
              % the new energy
              e2 = lnZ_q2 + lnZ_r2 -lnZ_s2;
              
              % check that the energy is defined and that the tilted moments are proper
              if ~all(p) || ~isfinite(e2)
                df=df*0.5;
                if df<0.1,
                  break
                else
                  if ismember(display,{'iter'})
                    fprintf('%d, e=%.6f, dm=%.4f, dV=%.4f, increasing damping to df=%g.\n',i1,e,tol_m(1),tol_m(2),df)
                  end
                  continue
                end
              end
              
              % accept the new state
              [nu_q,tau_q,mf,Vf,Sf,lnZ_q]=deal(nu_q2,tau_q2,mf2,Vf2,Sf2,lnZ_q2);
              [lnZ_r,lnZ_i,m_r,V_r,lnZ_s,nu_s,tau_s]=deal(lnZ_r2,lnZ_i2,m_r2,V_r2,lnZ_s2,nu_s2,tau_s2);
              
              % EP search direction (moment matching)
              [dnu_q,dtau_q]=ep_update_dir(mf,Vf,m_r,V_r,eta,up_mode,tolUpdate);
              
              % Check for convergence
              % the difference between the marginal moments
%               Vf=diag(Sf);
              tol_m=[abs(mf-m_r) abs(Vf-V_r)];
              
              % measure the convergence by the moment difference
              convergence=all(tol_m(:,1)<tolStop*abs(mf)) && all(tol_m(:,2)<tolStop*abs(Vf));
              
              % measure the convergence by the change of energy
              %convergence=abs(e2-e)<tolStop;
              
              tol_m=max(tol_m);
              e=e2;
              
              if ismember(display,{'iter'})
                fprintf('%d, e=%.6f, dm=%.4f, dV=%.4f, df=%g.\n',i1,e,tol_m(1),tol_m(2),df)
              end
              
              if convergence
                if ismember(display,{'final','iter'})
                  fprintf('Convergence with parallel EP, iter %d, e=%.6f, dm=%.4f, dV=%.4f, df=%g.\n',i1,e,tol_m(1),tol_m(2),df)
                end
                break
              end
            end
          end % end of initial rounds of parallel EP
          
          if isfinite(e) && ~convergence
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if no convergence with the parallel EP
            % start double-loop iterations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            up_mode=gp.latent_opt.up_mode; % update mode in double-loop iterations
                                           %up_mode='ep'; % choose the moment matching
                                           %up_mode='grad'; % choose the gradients
            df_lim=gp.latent_opt.df_lim; % step size limit (1 suitable for ep updates)
            
            tol_e=inf;  % the energy difference for measuring convergence (tol_e < tolStop)
            ninner=0;   % counter for the inner loop iterations
            df=df0;     % initial step size (damping factor)
            
            % the intial gradient in the search direction
            g = sum( (mf -m_r).*dnu_q ) +0.5*sum( (V_r +m_r.^2 -Vf -mf.^2).*dtau_q );
            
            sdir_reset=false;
            rec_sadj=[0 e g]; % record for step size adjustment
            for i1=1:maxiter
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % calculate a new proposal state
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
              % Limit the step size separately for each site so that the cavity variances
              % do not exceed the upper limit (this will change the search direction)
              % this should not happen after step size adjustment
              ii1=tau_s-eta.*(tau_q+df*dtau_q)<tauc_min;
              if any(ii1)
                %ii1=dtau_q>0; df1=min( ( (tau_s(ii1)-tauc_min(ii1))./eta(ii1)-tau_q(ii1) )./dtau_q(ii1)/df ,1);
                df1=( (tau_s(ii1)-tauc_min(ii1))./eta(ii1) -tau_q(ii1) )./dtau_q(ii1)/df;
                
                dnu_q(ii1)=dnu_q(ii1).*df1;
                dtau_q(ii1)=dtau_q(ii1).*df1;
                
                % the intial gradient in the search direction
                g = sum( (mf -m_r).*dnu_q ) +0.5*sum( (V_r +m_r.^2 -Vf -mf.^2).*dtau_q );
                
                % re-init the step size adjustment record
                rec_sadj=[0 e g];
              end
              % proposal
              nu_q2=nu_q+df*dnu_q;
              tau_q2=tau_q+df*dtau_q;
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % energy for the proposal state
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
              % update the q-distribution
%               [mf2,Sf2,lnZ_q2,L1,L2]=evaluate_q(nu_q2,tau_q2,K,display,K_uu, K_fu, Kv_ff, Qv_ff);
              switch gp.type
                case 'FULL'
                  [mf2,Sf2,lnZ_q2,L1,L2]=evaluate_q(nu_q2,tau_q2,K,display);
                  Vf2 = diag(Sf2);
                case 'FIC'
                  [mf2,Vf2,lnZ_q2,L1,L2]=evaluate_q2(nu_q2,tau_q2,Luu, K_fu, Kv_ff, Qv_ff, display);
                otherwise
                  error('Robust-EP not implemented for this type of GP!');
              end
              
              % check cavity
              pcavity=all( (1./Vf2-eta.*tau_q2 )>=tauc_min);
              
              g2=NaN;
              if isempty(L2)
                % the q-distribution not defined (the posterior covariance
                % not positive definite)
                e2=inf;
              elseif pcavity
                % the tilted distribution
                [lnZ_r2,lnZ_i2,m_r2,V_r2]=evaluate_r(nu_q2,tau_q2,eta,fh_tm,nu_s,tau_s,display);
                
                % the new energy
                e2 = lnZ_q2 + lnZ_r2 -lnZ_s;
                
                % gradients in the search direction
                g2 = sum( (mf2 -m_r2).*dnu_q ) +0.5*sum( (V_r2 +m_r2.^2 -Vf2 -mf2.^2).*dtau_q );
                
                if ismember(display,{'iter'})
                  % ratio of the gradients
                  fprintf('dg=%6.3f, ',min(abs(g2)/abs(g),99))
                end
              end
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % check if the energy decreases
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if ~isfinite(e2) || ( pcavity && g2>10*abs(g) )
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % ill-conditioned q-distribution or very large increase
                % in the gradient
                % => half the step size
                df=df*0.5;
                
                if ismember(display,{'iter'})
                  fprintf('decreasing step size, ')
                end
              elseif ~pcavity && ~pvis
                % The cavity distributions resulting from the proposal distribution
                % are not well defined, reset the site parameters by doing 
                % one parallel update with a zero initialization and continue 
                % with double loop iterations
                
                if ismember(display,{'iter'})
                  fprintf('re-init the posterior due to ill-conditioned cavity distributions, ')
                end
                
                % Do resetting only once
                pvis=1;
                
                up_mode='ep';
                nu_q=zeros(size(y));tau_q=zeros(size(y));
                mf=zeros(size(y));
                switch gp.type
                  case 'FULL'
                    Sf=K;Vf=diag(K);
                  case 'FIC'
                    Vf=Cv_ff;
                end
                nu_s=mf./Vf;
                tau_s=1./Vf;
%                 lnZ_s=0.5*sum( (-log(tau_s) +nu_s.^2 ./tau_s)./eta ); % minus 0.5*log(2*pi)./eta
                [lnZ_r,lnZ_i,m_r,V_r]=evaluate_r(nu_q,tau_q,eta,fh_tm,nu_s,tau_s,display);
%                 e = lnZ_q + lnZ_r -lnZ_s;
                [dnu_q,dtau_q]=ep_update_dir(mf,Vf,m_r,V_r,eta,up_mode,tolUpdate);
                %nu_q=dnu_q; tau_q=dtau_q;
                nu_q=0.9.*dnu_q; tau_q=0.9.*dtau_q;
                
                switch gp.type
                  case 'FULL'
                    [mf,Sf,lnZ_q]=evaluate_q(nu_q,tau_q,K,display);
                    Vf = diag(Sf);
                  case 'FIC'
                    [mf,Vf,lnZ_q]=evaluate_q2(nu_q,tau_q,Luu, K_fu, Kv_ff, Qv_ff, display);
                  otherwise
                    error('Robust-EP not implemented for this type of GP!');
                end
                nu_s=mf./Vf; tau_s=1./Vf;
                lnZ_s=0.5*sum( (-log(tau_s) +nu_s.^2 ./tau_s)./eta ); % minus 0.5*log(2*pi)./eta
                [lnZ_r,lnZ_i,m_r,V_r]=evaluate_r(nu_q,tau_q,eta,fh_tm,nu_s,tau_s,display);
                e = lnZ_q + lnZ_r -lnZ_s;
                [dnu_q,dtau_q]=ep_update_dir(mf,Vf,m_r,V_r,eta,up_mode,tolUpdate);
                
                df=0.8;
                
                g = sum( (mf -m_r).*dnu_q ) +0.5*sum( (V_r +m_r.^2 -Vf -mf.^2).*dtau_q );
                rec_sadj=[0 e g];
                
              elseif size(rec_sadj,1)<=1 && ( e2>e || abs(g2)>abs(g)*tolGrad )
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % no decrease in energy or the new gradient exceeds the
                % pre-defined limit
                % => adjust the step size
                
                if ismember(display,{'iter'})
                  fprintf('adjusting step size,  ')
                end
                
                % update the record for step size adjustment
                ii1=find(df>rec_sadj(:,1),1,'last');
                ii2=find(df<rec_sadj(:,1),1,'first');
                rec_sadj=[rec_sadj(1:ii1,:); df e2 g2; rec_sadj(ii2:end,:)];
                
                df_new=0;
                if size(rec_sadj,1)>1
                  if exist('csape','file')==2
                    if g2>0 
                      % adjust the step size with spline interpolation
                      pp=csape(rec_sadj(:,1)',[rec_sadj(1,3) rec_sadj(:,2)' rec_sadj(end,3)],[1 1]);
                      [tmp,df_new]=fnmin(pp,[0 df]);
                      
                    elseif isfinite(g2)
                      % extrapolate with Hessian end-conditions
                      H=(rec_sadj(end,3)-rec_sadj(end-1,3))/(rec_sadj(end,1)-rec_sadj(end-1,1));
                      pp=csape(rec_sadj(:,1)',[rec_sadj(1,3) rec_sadj(:,2)' H],[1 2]);
                      % extrapolate at most by 100% at a time
                      [tmp,df_new]=fnmin(pp,[df df*1.5]);
                    end
                  else
                    % if curvefit toolbox does not exist, use a simple Hessian
                    % approximation
                    [tmp,ind]=sort(rec_sadj(:,2),'ascend');
                    ind=ind(1:2);
                    
                    H=(rec_sadj(ind(1),3)-rec_sadj(ind(2),3))/(rec_sadj(ind(1),1)-rec_sadj(ind(2),1));
                    df_new=rec_sadj(ind(1),1) -rec_sadj(ind(1),3)/H;
                    if g2>0
                      % interpolate
                      df_new=max(min(df_new,df),0);
                    else
                      % extrapolate at most 100%
                      df_new=max(min(df_new,1.5*df),df);
                    end
                  end
                  df_new=min(df_new,df_lim);
                end
                
                if df_new==0
                  % the spline approxmation fails or no record of the previous gradients
                  if g2>0
                    df=df*0.9; % too long step since the gradient is positive
                  else
                    df=df*1.1; % too short step since the gradient is negative
                  end
                else
                  df=df_new;
                end
                % prevent too small cavity-variances after the step-size adjustment
                ii1=dtau_q>0;
                if any(ii1)
                  df_max=min( ( (tau_s(ii1)-tauc_min(ii1)-1e-8)./eta(ii1) -tau_q(ii1) )./dtau_q(ii1) );
                  df=min(df,df_max);
                end
                
              elseif e2>e+tolInner || (abs(g2)>abs(g)*tolGrad && strcmp(up_mode,'ep'))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % No decrease in energy despite the step size adjustments.
                % In some difficult cases the EP search direction may not
                % result in decrease of the energy or the gradient
                % despite of the step size adjustment. One reason for this
                % may be the parallel EP search direction
                % => try the negative gradient as the search direction
                %
                % or if the problem persists
                % => try resetting the search direction
                
                if abs(g2)>abs(g)*tolGrad && strcmp(up_mode,'ep')
                  % try switching to gradient based updates
                  up_mode='grad';
                  df_lim=1e3;
                  df=0.1;
                  if ismember(display,{'iter'})
                    fprintf('switch to gradient updates, ')
                  end
                elseif ~sdir_reset
                  if ismember(display,{'iter'})
                    fprintf('reset the search direction, ')
                  end
                  sdir_reset=true;
                elseif g2<0 && abs(g2)<abs(g) && e2>e
                  if ismember(display,{'final','iter'})
                    fprintf('Unable to continue: gradients of the inner-loop objective are inconsistent\n')
                  end
                  break;
                else
                  df=df*0.1;
                end
                
                % the new search direction
                [dnu_q,dtau_q]=ep_update_dir(mf,Vf,m_r,V_r,eta,up_mode,tolUpdate);
                
                % the initial gradient in the search direction
                g = sum( (mf -m_r).*dnu_q ) +0.5*sum( (V_r +m_r.^2 -Vf -mf.^2).*dtau_q );
                
                % re-init the step size adjustment record
                rec_sadj=[0 e g];
              else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % decrease of energy => accept the new state
                
                dInner=abs(e-e2); % the inner loop energy change
                
                % accept the new site parameters (nu_q,tau_q)
                [mf,Vf,Sf,nu_q,tau_q,lnZ_q]=deal(mf2,Vf2,Sf2,nu_q2,tau_q2,lnZ_q2);
                
                % accept also the new tilted distributions
                [lnZ_r,lnZ_i,m_r,V_r,e]=deal(lnZ_r2,lnZ_i2,m_r2,V_r2,e2);
                
                % check that the new cavity variances are positive and not too large
                tau_s2=1./Vf;
                pcavity=all( (tau_s2-eta.*tau_q )>=tauc_min);
                supdate=false;
                if pcavity && (dInner<tolInner || ninner>=max_ninner)
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % try to update the surrogate distribution on the condition that
                  % - the cavity variances are positive and not too large
                  % - the new tilted moments are proper
                  % - sufficient tolerance or the maximum number of inner
                  %   loop updates is exceeded
                  
                  % update the surrogate distribution
                  nu_s2=mf.*tau_s2;
                  lnZ_s2=0.5*sum( (-log(tau_s2) +nu_s2.^2 ./tau_s2)./eta );
                  
                  % update the tilted distribution
                  [lnZ_r2,lnZ_i2,m_r2,V_r2]=evaluate_r(nu_q,tau_q,eta,fh_tm,nu_s2,tau_s2,display);
                  
                  % evaluate the new energy
                  e2 = lnZ_q + lnZ_r2 -lnZ_s2;
                  
                  if isfinite(e2)
                    % a successful surrogate update
                    supdate=true;
                    ninner=0; % reset the inner loop iteration counter
                    
                    % update the convergence criteria
                    tol_e=abs(e2-e);
                    
                    % accept the new state
                    [lnZ_r,lnZ_i,m_r,V_r,lnZ_s,nu_s,tau_s,e]=deal(lnZ_r2,lnZ_i2,m_r2,V_r2,lnZ_s2,nu_s2,tau_s2,e2);
                    
                    if ismember(display,{'iter'})
                      fprintf('surrogate update,     ')
                    end
                  else
                    % Improper tilted moments even though the cavity variances are
                    % positive. This is an indication of numerically unstable
                    % tilted moment integrations but fractional updates usually help
                    % => try switching to fractional updates
                    pcavity=false;
                    
                    if ismember(display,{'iter'})
                      fprintf('surrogate update failed, ')
                    end
                  end
                end
                
                if all(eta==eta1) && ~pcavity && (dInner<tolInner || ninner>=max_ninner)
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % If the inner loop moments (within tolerance) are matched
                  % but the new cavity variances are negative or the tilted moment
                  % integrations fail after the surrogate update
                  % => switch to fractional EP.
                  %
                  % This is a rare situation and most likely the
                  % hyperparameters are such that the approximating family
                  % is not flexible enough, i.e., the hyperparameters are
                  % unsuitable for the data.
                  %
                  % One can also try to reduce the lower limit for the
                  % cavity precisions tauc_min=1./(Vc_lim*diag(K)), i.e.
                  % increase the maximum cavity variance Vc_lim.
                  
                  % try switching to fractional updates
                  eta=repmat(eta2,n,1);
                  
                  % correct the surrogate normalization accordingly
                  % the surrogate distribution is not updated
                  lnZ_s2=0.5*sum( (-log(tau_s) +nu_s.^2 ./tau_s)./eta );
                  
                  % update the tilted distribution
                  [lnZ_r2,lnZ_i2,m_r2,V_r2]=evaluate_r(nu_q,tau_q,eta,fh_tm,nu_s,tau_s,display);
                  
                  % evaluate the new energy
                  e2 = lnZ_q + lnZ_r2 -lnZ_s2;
                  
                  if isfinite(e2)
                    % successful switch to fractional energy
                    supdate=true;
                    pcavity=true;
                    ninner=0; % reset the inner loop iteration counter
                    
                    % accept the new state
                    [lnZ_r,lnZ_i,m_r,V_r,lnZ_s,e]=deal(lnZ_r2,lnZ_i2,m_r2,V_r2,lnZ_s2,e2);
                    
                    % start with ep search direction
                    up_mode='ep';
                    df_lim=0.9;
                    df=0.1;
                    if ismember(display,{'iter'})
                      fprintf('switching to fractional EP, ')
                    end
                  else
                    % Improper tilted moments even with fractional updates
                    % This is very unlikely to happen because decreasing the
                    % fraction parameter (eta2<eta1) stabilizes the
                    % tilted moment integrations
                    
                    % revert back to the previous fraction parameter
                    eta=repmat(eta1,n,1);
                    
                    if ismember(display,{'final','iter'})
                      fprintf('Unable to switch to the fractional EP, check that eta2<eta1\n')
                    end
                    break;
                  end
                end
                             
                if all(eta==eta2) && ~pcavity && (dInner<tolInner || ninner>=10)
                  % Surrogate updates do not result into positive cavity variances
                  % even with fractional updates with eta2 => terminate iterations
                  if ismember(display,{'final','iter'})
                    fprintf('surrogate update failed with fractional updates, try decreasing eta2\n')
                  end
                  break
                end
                
                if ~supdate
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % no successful surrogate update, no sufficient tolerance,
                  % or the maximum number of inner loop updates is not yet exceeded
                  % => continue with the same surrogate distribution
                  
                  ninner=ninner+1; % increase inner loop iteration counter
                  if ismember(display,{'iter'})
                    fprintf('inner-loop update,    ')
                  end
                end
                
                % the new search direction
                [dnu_q,dtau_q]=ep_update_dir(mf,Vf,m_r,V_r,eta,up_mode,tolUpdate);
                
                % the initial gradient in the search direction
                g = sum( (mf -m_r).*dnu_q ) +0.5*sum( (V_r +m_r.^2 -Vf -mf.^2).*dtau_q );
                
                % re-init step size adjustment record
                rec_sadj=[0 e g];
              end
              
              if ismember(display,{'iter'})
                % maximum difference of the marginal moments
                tol_m=[max(abs(mf-m_r)) max(abs(Vf-V_r))];
                fprintf('%d, e=%.6f, dm=%.4f, dV=%.4f, df=%6f, eta=%.2f\n',i1,e,tol_m(1),tol_m(2),df,eta(1))
              end
              
              %%%%%%%%%%%%%%%%%%%%%%%
              % check for convergence
              convergence = tol_e<=tolStop;
              if convergence
                if ismember(display,{'final','iter'})
                  % maximum difference of the marginal moments
                  tol_m=[max(abs(mf-m_r)) max(abs(Vf-V_r))];
                  fprintf('Convergence, iter %d, e=%.6f, dm=%.4f, dV=%.4f, df=%6f, eta=%.2f\n',i1,e,tol_m(1),tol_m(2),df,eta(1))
                end
                break
              end
            end % end of the double-loop updates
          end
          
          % the current energy is not finite or no convergence
          if ~isfinite(e)
            fprintf('GPEP_E: Initial energy not defined, check the hyperparameters\n')
          elseif ~convergence
            fprintf('GPEP_E: No convergence, %d iter, e=%.6f, dm=%.4f, dV=%.4f, df=%6f, eta=%.2f\n',i1,e,tol_m(1),tol_m(2),df,eta(1))
            fprintf('GPEP_E: Check the hyperparameters, increase maxiter and/or max_ninner, or decrease tolInner\n')
          end
          edata=-e; % the data contribution to the marginal posterior density
          
          % =====================================================================================
          % Evaluate the prior contribution to the error from covariance functions and likelihood
          % =====================================================================================
          
          % Evaluate the prior contribution to the error from covariance functions
          eprior = 0;
          for i=1:ncf
            gpcf = gp.cf{i};
            eprior = eprior - gpcf.fh.lp(gpcf);
            %         eprior = eprior - feval(gpcf.fh.lp, gpcf, x, y);
          end
          
          % Evaluate the prior contribution to the error from likelihood functions
          if isfield(gp, 'lik') && isfield(gp.lik, 'p')
            likelih = gp.lik;
            eprior = eprior - likelih.fh.lp(likelih);
          end
          
          % the total energy
          e = edata + eprior;
          
          sigm2vec_i = 1./(tau_s-eta.*tau_q);     % vector of cavity variances
          muvec_i = (nu_s-eta.*nu_q).*sigm2vec_i; % vector of cavity means
          logZ_i = lnZ_i; % vector of tilted normalization factors

          
          % check that the posterior covariance is positive definite and
          % calculate its Cholesky decomposition
          switch gp.type
            case 'FULL'
              [L, notpositivedefinite] = chol(Sf);
              b = [];
              La2 = [];
              if notpositivedefinite || ~isfinite(e)
                [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite();
                return
              end
            case 'FIC'
              La2 = Luu;
              L = L2;
              b = Kv_ff - Qv_ff;
          end

              
          nutilde = nu_q;
          tautilde = tau_q;
          
          % store values to the cache
          ch.w = w;
          ch.e = e;
          ch.edata = edata;
          ch.eprior = eprior;
          ch.L = L;
          ch.nu_q = nu_q;
          ch.tau_q = tau_q;
          ch.La2 = La2;
          ch.b = b;
          ch.eta = eta;
          ch.logZ_i = logZ_i;
          ch.sigm2vec_i = sigm2vec_i;
          ch.muvec_i = muvec_i;
          ch.datahash = datahash;
          
        end
      otherwise
        error('Unknown optim method!');
    end
  end

  function [e, edata, eprior, tautilde, nutilde, L, La2, b, muvec_i, sigm2vec_i, logZ_i, eta, ch] = set_output_for_notpositivedefinite()
    % Instead of stopping to chol error, return NaN
    e = NaN;
    edata = NaN;
    eprior = NaN;
    tautilde = NaN;
    nutilde = NaN;
    L = NaN;
    La2 = NaN;
    b = NaN;
    muvec_i = NaN;
    sigm2vec_i = NaN;
    logZ_i = NaN;
    datahash = NaN;
    eta = NaN;
    w = NaN;
    ch.e = e;
    ch.edata = edata;
    ch.eprior = eprior;
    ch.tautilde = tautilde;
    ch.nutilde = nutilde;
    ch.L = L;
    ch.La2 = La2;
    ch.b = b;
    ch.muvec_i = muvec_i;
    ch.sigm2vec_i = sigm2vec_i;
    ch.logZ_i = logZ_i;
    ch.eta = eta;
    ch.datahash=datahash;
    ch.w = NaN;
  end
    
end

function [m_q,S_q,lnZ_q,L1,L2]=evaluate_q(nu_q,tau_q,K,display)

% function for determining the parameters of the q-distribution
% when site variances tau_q may be negative
%
% q(f) = N(f|0,K)*exp( -0.5*f'*diag(tau_q)*f + nu_q'*f )/Z_q = N(f|m_q,S_q)
%
% S_q = inv(inv(K)+diag(tau_q))
% m_q = S_q*nu_q;
%
% det(eye(n)+K*diag(tau_q))) = det(L1)^2 * det(L2)^2
% where L1 and L2 are upper triangular
%
% see Expectation consistent approximate inference (Opper & Winther, 2005)

n=length(nu_q);
ii1=find(tau_q>0); n1=length(ii1); W1=sqrt(tau_q(ii1));
ii2=find(tau_q<0); n2=length(ii2); W2=sqrt(abs(tau_q(ii2)));

L=zeros(n);
S_q=K;
if ~isempty(ii1)
  % Cholesky decomposition for the positive sites
  L1=(W1*W1').*K(ii1,ii1);
  L1(1:n1+1:end)=L1(1:n1+1:end)+1;
  L1=chol(L1);
  
  L(:,ii1) = bsxfun(@times,K(:,ii1),W1')/L1;
  
  S_q=S_q-L(:,ii1)*L(:,ii1)';
else
  L1=1;
end

if ~isempty(ii2)
  % Cholesky decomposition for the negative sites
  V=bsxfun(@times,K(ii2,ii1),W1')/L1;
  L2=(W2*W2').*(V*V'-K(ii2,ii2));
  L2(1:n2+1:end)=L2(1:n2+1:end)+1;
  
  [L2,pd]=chol(L2);
  if pd==0
    L(:,ii2)=bsxfun(@times,K(:,ii2),W2')/L2 -L(:,ii1)*(bsxfun(@times,V,W2)'/L2);
    S_q=S_q+L(:,ii2)*L(:,ii2)';
  else
    L2=[];
    if ismember(display,{'iter'})
      fprintf('Negative definite q-distribution.\n')
    end
  end
  
else
  L2=1;
end
%V_q=diag(S_q);
m_q=S_q*nu_q;

% log normalization
lnZ_q = -sum(log(diag(L1))) -sum(log(diag(L2))) +0.5*sum(m_q.*nu_q);

end

function [m_q,S_q,lnZ_q,L1,L2]=evaluate_q2(nu_q,tau_q,LK_uu, K_fu, Kv_ff, Qv_ff, display)

% function for determining the parameters of the q-distribution
% when site variances tau_q may be negative
%
% q(f) = N(f|0,K)*exp( -0.5*f'*diag(tau_q)*f + nu_q'*f )/Z_q = N(f|m_q,S_q)
%
% S_q = inv(inv(K)+diag(tau_q)) where K is sparse approximation for prior
%       covariance
% m_q = S_q*nu_q;
%
% det(eye(n)+K*diag(tau_q))) = det(L1)^2 * det(L2)^2
% where L1 and L2 are upper triangular
%
% see Expectation consistent approximate inference (Opper & Winther, 2005)

n=length(nu_q);

S_q = Kv_ff;
m_q = nu_q;
D = Kv_ff - Qv_ff;
L1 = sqrt(1 + D.*tau_q);
L = [];
if any(~isreal(L1))
  if ismember(display,{'iter'})
    fprintf('Negative definite q-distribution.\n')
  end
else
  U = K_fu;
  WDtilde = tau_q./(1+tau_q.*D);
  
  % Evaluate diagonal of S_q
  
  ii1=find(WDtilde>0); n1=length(ii1); W1=sqrt(WDtilde(ii1)); % WS^-1
  ii2=find(WDtilde<0); n2=length(ii2); W2=sqrt(abs(WDtilde(ii2))); % WS^-1
  if ~isempty(ii2) || ~isempty(ii1)
    if ~isempty(ii1)
      UWS(:,ii1) = bsxfun(@times, U(ii1,:)', W1');
    end
    
    if ~isempty(ii2)
      UWS(:,ii2) = bsxfun(@times, U(ii2,:)', W2');
    end
    [L, p] = chol(LK_uu*LK_uu' + UWS(:,ii1)*UWS(:,ii1)' - UWS(:,ii2)*UWS(:,ii2)', 'lower');
    if p~=0
      L=[];
      if ismember(display,{'iter'})
        fprintf('Negative definite q-distribution.\n')
      end
    else
      
      S = 1 + D.*tau_q;
%               S_q = diag(D./S) + diag(1./S)*U*inv(L*L')*U'*diag(1./S);
      S_q = D./S + sum((bsxfun(@times, 1./S, U)/L').^2,2);
      m_q = D.*nu_q./S + (U*(L'\(L\(U'*(nu_q./S)))))./S;
    end
  else
  end
  %   end
  
end

% log normalization
L2 = L;
lnZ_q = -0.5*sum(log(L1.^2)) - sum(log(diag(L))) + sum(log(diag(LK_uu))) +0.5*sum(m_q.*nu_q);

end

function [lnZ_r,lnZ_i,m_r,V_r,p]=evaluate_r(nu_q,tau_q,eta,fh_tm,nu_s,tau_s,display)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for determining the parameters of the r-distribution 
% (the product of the tilted distributions)
%
% r(f) = exp(-lnZ_r) * prod_i p(y(i)|f(i)) * exp( -0.5*f(i)^2 tau_r(i) + nu_r(i)*f(i) )
%      ~ prod_i N(f(i)|m_r(i),V_r(i))
%
% tau_r = tau_s - tau_q
% nu_r = nu_s - nu_q
%
% lnZ_i(i) = log int p(y(i)|f(i)) * N(f(i)|nu_r(i)/tau_r(i),1/tau_r(i)) df(i)
%
% see Expectation consistent approximate inference (Opper & Winther, 2005)

n=length(nu_q);
[lnZ_i,m_r,V_r,nu_r,tau_r]=deal(zeros(n,1));
p=false(n,1);
for si=1:n
  % cavity distribution
  tau_r_si=tau_s(si)-eta(si)*tau_q(si);
  if tau_r_si<=0
    %     if ismember(display,{'iter'})
    %       %fprintf('Negative cavity precision at site %d\n',si)
    %     end
    continue
  end
  nu_r_si=nu_s(si)-eta(si)*nu_q(si);
  
  % tilted moments
  [lnZ_si,m_r_si,V_r_si] = fh_tm(si, nu_r_si/tau_r_si, 1/tau_r_si, eta(si));
  
  if ~isfinite(lnZ_si) || V_r_si<=0
    %     if ismember(display,{'iter'})
    %       fprintf('Improper normalization or tilted variance at site %d\n',si)
    %     end
    continue
  end
  
  % store the new parameters
  [nu_r(si),tau_r(si),lnZ_i(si),m_r(si),V_r(si)]=deal(nu_r_si,tau_r_si,lnZ_si,m_r_si,V_r_si);
  
  p(si)=true;
end

lnZ_r=sum(lnZ_i./eta) +0.5*sum((-log(tau_r) +nu_r.^2 ./tau_r)./eta);
end

function [dnu_q,dtau_q]=ep_update_dir(m_q,V_q,m_r,V_r,eta,up_mode,tolUpdate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update direction for double-loop EP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% V_q=diag(S_q);
switch up_mode
  case 'ep'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % site updates by moment matching
    
    [dnu_q,dtau_q]=deal(zeros(size(m_q)));
    
    %ind_up=V_r>0 & max(abs(V_r-V_q),abs(m_r-m_q))>tolUpdate;
    ind_up=V_r>0 & (abs(V_r-V_q) > tolUpdate*abs(V_q) | abs(m_r-m_q) > tolUpdate*abs(m_q));
    
    dnu_q(ind_up) = ( m_r(ind_up)./V_r(ind_up) - m_q(ind_up)./V_q(ind_up) ) ./ eta(ind_up);
    dtau_q(ind_up) = ( 1./V_r(ind_up) - 1./V_q(ind_up) )./ eta(ind_up);
    
  case 'grad'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % gradient descend
    % Not used at the moment!
    
    % evaluate the gradients wrt nu_q and tau_q
    gnu_q = m_q - m_r;
    gtau_q = 0.5*(V_r + m_r.^2 - V_q - m_q.^2);
    
    % the search direction
    dnu_q=-gnu_q;
    dtau_q=-gtau_q;  
end

end
