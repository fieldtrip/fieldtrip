function [g, gdata, gprior] = gpep_g(w, gp, x, y, varargin)
%GPEP_G  Evaluate gradient of EP's marginal log posterior estimate (GPEP_E)
%
%  Description
%    G = GPEP_G(W, GP, X, Y, OPTIONS) takes a full GP parameter
%    vector W, structure GP, a matrix X of input vectors and a
%    matrix Y of target vectors, and evaluates the gradient G of
%    EP's marginal log posterior estimate (gpep_e). Each row of X
%    corresponds to one input vector and each row of Y corresponds
%    to one target vector.
%
%    [G, GDATA, GPRIOR] = GPEP_G(GP, X, Y, OPTIONS) also returns
%    separately the data and prior contributions to the gradient.
%    
%    OPTIONS is optional parameter-value pair
%      z - optional observed quantity in triplet (x_i,y_i,z_i)
%          Some likelihoods may use this. For example, in case of
%          Poisson likelihood we have z_i=E_i, that is, expected
%          value for ith case.
%
%  See also
%    GP_SET, GP_G, GPEP_E, GPEP_PRED

% Copyright (c) 2007-2010  Jarno Vanhatalo
% Copyright (c) 2010       Heikki Peura, Aki Vehtari
  
% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.
  
  ip=inputParser;
  ip.FunctionName = 'GPEP_G';
  ip.addRequired('w', @(x) isvector(x) && isreal(x) && all(isfinite(x)));
  ip.addRequired('gp',@isstruct);
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addOptional('method', '1', @(x) ismember(x,{'1','2'}))
  ip.parse(w, gp, x, y, varargin{:});
  z=ip.Results.z;
  method = ip.Results.method;

  
  gp=gp_unpak(gp, w);       % unpak the parameters
  ncf = length(gp.cf);
  n=size(x,1);
  
  g = [];
  gdata = [];
  gprior = [];
  
  if isfield(gp,'savememory') && gp.savememory
    savememory=1;
  else
    savememory=0;
  end
  
  % First Evaluate the data contribution to the error    
  switch gp.type
    % ============================================================
    % FULL
    % ============================================================
    case 'FULL'   % A full GP
                  % Calculate covariance matrix and the site parameters
      [K, C] = gp_trcov(gp,x);
      
      if issparse(C)          
        % If compact support covariance functions are used 
        % the covariance matrix will be sparse
        [e, edata, eprior, tautilde, nutilde, LD, tmp, tmp, mu_i, sigm2_i, Z_i, eta] = gpep_e(w, gp, x, y, 'z', z);
        Stildesqroot = sparse(1:n,1:n,sqrt(tautilde),n,n);
        
        b = nutilde - Stildesqroot*ldlsolve(LD,Stildesqroot*(C*nutilde));
        % evaluate the sparse inverse
        invC = spinv(LD,1);       
        invC = Stildesqroot*invC*Stildesqroot;
      else
        [e, edata, eprior, tautilde, nutilde, L, tmp, tmp, mu_i, sigm2_i, Z_i, eta] = gpep_e(w, gp, x, y, 'z', z);
        
        if all(tautilde > 0) && ~isequal(gp.latent_opt.optim_method, 'robust-EP')
          % This is the usual case where likelihood is log concave
          % for example, Poisson and probit
          % Stildesqroot=diag(sqrt(tautilde));
          % logZep; nutilde; tautilde;
          % b=nutilde-Stildesqroot*(L'\(L\(Stildesqroot*(C*nutilde))));
          % invC = Stildesqroot*(L'\(L\Stildesqroot));
          Stildesqroot=sqrt(tautilde);
          temp=L\diag(Stildesqroot);
          invC = temp'*temp;
          b=nutilde-Stildesqroot.*(L'\(L\(Stildesqroot.*(C*nutilde))));
        elseif isequal(gp.latent_opt.optim_method, 'robust-EP')
          
          A=bsxfun(@times,tautilde,L');   % Sf = L'*L;
          b=nutilde-A*(L*nutilde);        % (eye(n)-diag(tautilde)*Sf)\nutilde
          A=-A*A';                         % = -diag(tautilde)*Sf*diag(tautilde)
          A(1:n+1:end)=A(1:n+1:end)+tautilde'; % diag(tautilde)-diag(tautilde)*Sf*diag(tautilde)
          invC = A;

        else                         
          % We might end up here if the likelihood is not log concace
          % For example Student-t likelihood.
          % NOTE! This does not work reliably yet
          S = diag(tautilde);
          b = nutilde - tautilde.*(L'*L*(nutilde));
          invC = S*L';
          invC = S - invC*invC';
        end           
        
      end

      if ~all(isfinite(e));
        % instead of stopping to error, return NaN
        g=NaN;
        gdata = NaN;
        gprior = NaN;
        return;
      end
      
      % =================================================================
      % Gradient with respect to covariance function parameters
      if ~isempty(strfind(gp.infer_params, 'covariance'))        
        % Evaluate the gradients from covariance functions
        for i=1:ncf
          i1=0;
          if ~isempty(gprior)
            i1 = length(gprior);
          end
          gpcf = gp.cf{i};
          if savememory
            np=gpcf.fh.cfg(gpcf,[],[],[],0);
          else
            DKffc = gpcf.fh.cfg(gpcf, x);
            np=length(DKffc);
          end
          gprior_cf = -gpcf.fh.lpg(gpcf);
          
          if ~isfield(gp,'meanf')
            for i2 = 1:np
              if savememory
                DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
              else
                DKff=DKffc{i2};
              end
              i1 = i1+1;                
              Bdl = b'*(DKff*b);              
              Cdl = sum(sum(invC.*DKff)); % help arguments for lengthScale
              gdata(i1)=0.5.*(Cdl - Bdl);
              gprior(i1) = gprior_cf(i2);
            end
          else
            i1=0;
            Stildesqroot=diag(sqrt(tautilde));
            invKs=eye(size(C))-Stildesqroot*(L'\(L\(Stildesqroot*C)));
            [dMNM trA]=mean_gf(gp,x,C,invKs,DKff,Stildesqroot,nutilde,'EP');
            for i2 = 1:np
              i1=i1+1;
              if savememory
                DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
              else
                DKff=DKffc{i2};
              end
              trK=sum(sum(invC.*DKff));
              gdata(i2)=0.5*(-1*dMNM + trK + trA{i2});
              gprior(i1) = gprior_cf(i2);
            end
          end
          
          % Set the gradients of hyperparameter
          if length(gprior_cf) > np
            for i2=np+1:length(gprior_cf)
              i1 = i1+1;
              gdata(i1) = 0;
              gprior(i1) = gprior_cf(i2);
            end
          end
        end
      end
      
      % =================================================================
      % Gradient with respect to likelihood function parameters
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh, 'siteDeriv')
%         [Ef, Varf] = gpep_pred(gp, x, y, x, 'z', z);
        
%         sigm2_i = (Varf.^-1 - tautilde).^-1;
%         mu_i = sigm2_i.*(Ef./Varf - nutilde);
        
        gdata_lik = 0;
        lik = gp.lik;
        for k1 = 1:length(y)
          if isempty(eta)
            gdata_lik = gdata_lik - lik.fh.siteDeriv(lik, y, k1, sigm2_i(k1), mu_i(k1), z);
          else
            gdata_lik = gdata_lik - lik.fh.siteDeriv2(lik, y, k1, sigm2_i(k1), mu_i(k1), z, eta(k1), Z_i(k1));
          end
        end
        % evaluate prior contribution for the gradient
        if isfield(gp.lik, 'p')
          g_logPrior = -lik.fh.lpg(lik);
        else
          g_logPrior = zeros(size(gdata_lik));
        end
        % set the gradients into vectors that will be returned
        gdata = [gdata gdata_lik];
        gprior = [gprior g_logPrior];
        i1 = length(gdata);
      end
      
    case {'FIC'}
      % ============================================================
      % FIC
      % ============================================================
      g_ind = zeros(1,numel(gp.X_u));
      gdata_ind = zeros(1,numel(gp.X_u));
      gprior_ind = zeros(1,numel(gp.X_u));
      
      u = gp.X_u;
      DKuu_u = 0;
      DKuf_u = 0;

      [e, edata, eprior, tautilde, nutilde, L, La, b, mu_i, sigm2_i, Z_i, eta] = gpep_e(w, gp, x, y, 'z', z);

      K_fu = gp_cov(gp, x, u);         % f x u
      K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu        
      iKuuKuf = K_uu\K_fu';
      
      LL = sum(L.*L,2);
      
      if isequal(gp.latent_opt.optim_method, 'robust-EP')
        % Help parameters for Robust-EP
        S = 1+tautilde.*b;
        WiS = tautilde./S;
      end
      
      if ~all(isfinite(e));
        % instead of stopping to error, return NaN
        g=NaN;
        gdata = NaN;
        gprior = NaN;
        return;
      end
      
      % =================================================================
      % Gradient with respect to covariance function parameters
      if ~isempty(strfind(gp.infer_params, 'covariance'))
        for i=1:ncf            
          i1=0;
          if ~isempty(gprior)
            i1 = length(gprior);
          end
          
          gpcf = gp.cf{i};
          % Get the gradients of the covariance matrices 
          % and gprior from gpcf_* structures
          gpcf = gp.cf{i};
          if savememory
            % If savememory option is used, just get the number of
            % hyperparameters and calculate gradients later
            np=gpcf.fh.cfg(gpcf,[],[],[],0);
          else
            DKffc = gpcf.fh.cfg(gpcf, x, [], 1);
            DKuuc = gpcf.fh.cfg(gpcf, u);
            DKufc = gpcf.fh.cfg(gpcf, u, x);
            np=length(DKuuc);
          end
          gprior_cf = -gpcf.fh.lpg(gpcf);
          
          if ~isequal(gp.latent_opt.optim_method, 'robust-EP')
            
            for i2 = 1:np
              i1 = i1+1;
              if savememory
                DKff=gpcf.fh.cfg(gpcf,x,[],1,i2);
                DKuu=gpcf.fh.cfg(gpcf,u,[],[],i2);
                DKuf=gpcf.fh.cfg(gpcf,u,x,[],i2);
              else
                DKff=DKffc{i2};
                DKuu=DKuuc{i2};
                DKuf=DKufc{i2};
              end
              KfuiKuuKuu = iKuuKuf'*DKuu;
              gdata(i1) = -0.5.*((2*b*DKuf'-(b*KfuiKuuKuu))*(iKuuKuf*b') + 2.*sum(sum(L'.*(L'*DKuf'*iKuuKuf))) - ...
                sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
              
              gdata(i1) = gdata(i1) - 0.5.*(b.*DKff')*b';
              gdata(i1) = gdata(i1) + 0.5.*(2.*b.*sum(DKuf'.*iKuuKuf',2)'*b'- b.*sum(KfuiKuuKuu.*iKuuKuf',2)'*b');
              gdata(i1) = gdata(i1) + 0.5.*(sum(DKff./La) - sum(LL.*DKff));
              gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf'.*iKuuKuf',2)) - sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
              
              gprior(i1) = gprior_cf(i2);
            end
            
          else
            % Robust-EP
            for i2 = 1:np
              i1 = i1+1;
              
              if savememory
                DKff=gpcf.fh.cfg(gpcf,x,[],1,i2);
                DKuu=gpcf.fh.cfg(gpcf,u,[],[],i2);
                DKuf=gpcf.fh.cfg(gpcf,u,x,[],i2);
              else
                DKff=DKffc{i2};
                DKuu=DKuuc{i2};
                DKuf=DKufc{i2};
              end
              % Evaluate derivative of log(det(I+Ktilde*W)) where Ktilde is
              % FIC sparse approximation of covariance with respect to
              % hyperparameters
              
              Dd = DKff{i2} - 2.*sum(DKuf.*iKuuKuf)' - sum((-iKuuKuf'*DKuu)'.*iKuuKuf)'; % d(diag(Kff - Qff)) / dth
              DS = Dd.*tautilde;
              gdata(i1) = -0.5.*sum(DS./S);
              DTtilde = DKuu + DKuf*bsxfun(@times, WiS, K_fu) - K_fu'*bsxfun(@times, WiS.*DS./S, K_fu) + ...
                K_fu'*bsxfun(@times, WiS, DKuf');
              gdata(i1) = gdata(i1) - 0.5.*sum(sum(inv(L).*(L\DTtilde)));
              gdata(i1) = gdata(i1) - 0.5.*sum(sum(-inv(La).*(La\DKuu))); 
              iSKfuiL = bsxfun(@times, 1./S, K_fu/L');
              
              % Evaluate derivative of quadratic term
              % nutilde'*sigma^-1*nutilde with respect to hyperparameters
            
              nud=(nutilde'*iSKfuiL)/L;
              nuDpcovnu = sum(nutilde.^2.*(-DS.*b./S.^2 + Dd./S)) + 2*(nutilde./S)'*(DKuf'*nud') - ...
                2*((nutilde.*DS./S)'*iSKfuiL)*(iSKfuiL'*nutilde) - nud*DTtilde*nud'; % nutilde'* d(sigma^-1)/dth *nutilde
              gdata(i1) = gdata(i1) + 0.5*nuDpcovnu; 
              gdata(i1) = -gdata(i1);
              gprior(i1) = gprior_cf(i2);
            end
            
          end
          
          % Set the gradients of hyperparameter
          if length(gprior_cf) > np
            for i2=np+1:length(gprior_cf)
              i1 = i1+1;
              gdata(i1) = 0;
              gprior(i1) = gprior_cf(i2);
            end
          end
        end
        
      end
      % =================================================================
      % Gradient with respect to a likelihood function parameters        
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh, 'siteDeriv')
%         [Ef, Varf] = gpep_pred(gp, x, y, x, 'tstind', 1:n, 'z', z);
%         sigm2_i = (Varf.^-1 - tautilde).^-1;
%         mu_i = sigm2_i.*(Ef./Varf - nutilde);
        
        gdata_lik = 0;
        lik = gp.lik;
        for k1 = 1:length(y)
          if ~isequal(gp.latent_opt.optim_method, 'robust-EP')
            gdata_lik = gdata_lik - lik.fh.siteDeriv(lik, y, k1, sigm2_i(k1), mu_i(k1), z);
          else
            gdata_lik = gdata_lik - lik.fh.siteDeriv2(lik, y, k1, sigm2_i(k1), mu_i(k1), z, eta(k1), Z_i(k1));
          end
        end

        % evaluate prior contribution for the gradient
        if isfield(gp.lik, 'p')
          g_logPrior = -lik.fh.lpg(lik);
        else
          g_logPrior = zeros(size(gdata_lik));
        end
        % set the gradients into vectors that will be returned
        gdata = [gdata gdata_lik];
        gprior = [gprior g_logPrior];
        i1 = length(gdata);
      end
      
      
      % =================================================================
      % Gradient with respect to inducing inputs
      
      if ~isempty(strfind(gp.infer_params, 'inducing'))
        if isfield(gp.p, 'X_u') && ~isempty(gp.p.X_u)
          m = size(gp.X_u,2);
          st=0;
          if ~isempty(gprior)
            st = length(gprior);
          end
          
          gdata(st+1:st+length(gp.X_u(:))) = 0;
          i1 = st+1;
          for i = 1:size(gp.X_u,1)
            if iscell(gp.p.X_u) % Own prior for each inducing input
              pr = gp.p.X_u{i};
              gprior(i1:i1+m) = pr.fh.lpg(gp.X_u(i,:), pr);
            else % One prior for all inducing inputs
              gprior(i1:i1+m-1) = gp.p.X_u.fh.lpg(gp.X_u(i,:), gp.p.X_u);
            end
            i1 = i1 + m;
          end
          
          for i=1:ncf
            i1=st;
            
            gpcf = gp.cf{i};
            if savememory
              % If savememory option is used, just get the number of
              % covariates in X and calculate gradients later
              np=gpcf.fh.ginput(gpcf,u,[],0);
            else
              DKuu = gpcf.fh.ginput(gpcf, u);
              DKuf = gpcf.fh.ginput(gpcf, u, x);
              np=1;
            end
            
            if ~isequal(gp.latent_opt.optim_method, 'robust-EP')
              for i3 = 1:np
                if savememory
                  DKuu=gpcf.fh.ginput(gpcf,u,[],i3);
                  DKuf=gpcf.fh.ginput(gpcf,u,x,i3);
                end
                for i2 = 1:length(DKuu)
                  i1 = i1+1;
                  
                  KfuiKuuKuu = iKuuKuf'*DKuu{i2};
                  
                  gdata(i1) = gdata(i1) - 0.5.*((2*b*DKuf{i2}'-(b*KfuiKuuKuu))*(iKuuKuf*b') + ...
                    2.*sum(sum(L'.*(L'*DKuf{i2}'*iKuuKuf))) - sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
                  gdata(i1) = gdata(i1) + 0.5.*(2.*b.*sum(DKuf{i2}'.*iKuuKuf',2)'*b'- b.*sum(KfuiKuuKuu.*iKuuKuf',2)'*b');
                  gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf{i2}'.*iKuuKuf',2)) - ...
                    sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
                end
              end
            else
              for i3 = 1:np
                if savememory
                  DKuu=gpcf.fh.ginput(gpcf,u,[],i3);
                  DKuf=gpcf.fh.ginput(gpcf,u,x,i3);
                end
                % Robust-EP
                for i2 = 1:length(DKuu)
                  i1 = i1+1;
                  Dd = -2.*sum(DKuf{i2}.*iKuuKuf)' - sum((-iKuuKuf'*DKuu{i2})'.*iKuuKuf)';
                  DS = Dd.*tautilde;
                  gdata(i1) = -0.5.*sum(DS./S);
                  DTtilde = DKuu{i2} + DKuf{i2}*bsxfun(@times, WiS, K_fu) - K_fu'*bsxfun(@times, WiS.*DS./S, K_fu) + ...
                    K_fu'*bsxfun(@times, WiS, DKuf{i2}');
                  gdata(i1) = gdata(i1) - 0.5.*sum(sum(inv(L).*(L\DTtilde)));
                  gdata(i1) = gdata(i1) - 0.5.*sum(sum(-La^-1.*(La\DKuu{i2})));
                  iSKfuiL = bsxfun(@times, 1./S, K_fu/L');
                  nud=(nutilde'*iSKfuiL)/L;
                  nuDpcovnu = sum(nutilde.^2.*(-DS.*b./S.^2 + Dd./S)) + 2*(nutilde./S)'*(DKuf{i2}'*nud') - 2*((nutilde.*DS./S)'*iSKfuiL)*(iSKfuiL'*nutilde) - nud*DTtilde*nud';
                  gdata(i1) = gdata(i1) + 0.5*nuDpcovnu;
                  gdata(i1) = -gdata(i1);
                end
              end
            end
          end
        end
      end
      

      
    case {'PIC' 'PIC_BLOCK'}
      % ============================================================
      % PIC
      % ============================================================
      g_ind = zeros(1,numel(gp.X_u));
      gdata_ind = zeros(1,numel(gp.X_u));
      gprior_ind = zeros(1,numel(gp.X_u));
      
      u = gp.X_u;
      ind = gp.tr_index;
      DKuu_u = 0;
      DKuf_u = 0;

      [e, edata, eprior, tautilde, nutilde, L, La, b] = gpep_e(w, gp, x, y, 'z', z);

      K_fu = gp_cov(gp, x, u);         % f x u
      K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu        
      iKuuKuf = K_uu\K_fu';
      
      % =================================================================
      % Gradient with respect to covariance function parameters
      if ~isempty(strfind(gp.infer_params, 'covariance'))
        
        % Evaluate the gradients from covariance functions
        for i=1:ncf            
          i1=0;
          if ~isempty(gprior)
            i1 = length(gprior);
          end
          
          % Get the gradients of the covariance matrices 
          % and gprior from gpcf_* structures
          gpcf = gp.cf{i};
          if savememory
            % If savememory option is used, just get the number of
            % hyperparameters and calculate gradients later
            np=gpcf.fh.cfg(gpcf,[],[],[],0);
          else
            DKuuc = gpcf.fh.cfg(gpcf, u);
            DKufc = gpcf.fh.cfg(gpcf, u, x);
            for kk = 1:length(ind)
              DKffc{kk} = gpcf.fh.cfg(gpcf, x(ind{kk},:));
            end
            np=length(DKuuc);
          end
          gprior_cf = -gpcf.fh.lpg(gpcf);
          
          for i2 = 1:np
            i1 = i1+1;
            if savememory
              DKuu=gpcf.fh.cfg(gpcf,u,[],[],i2);
              DKuf=gpcf.fh.cfg(gpcf,u,x,[],i2);
            else
              DKuu=DKuuc{i2};
              DKuf=DKufc{i2};
            end
            
            KfuiKuuKuu = iKuuKuf'*DKuu;
            %            H = (2*K_uf'- KfuiKuuKuu)*iKuuKuf;
            % Here we evaluate  gdata = -0.5.* (b*H*b' + trace(L*L'H)
            gdata(i1) = -0.5.*((2*b*DKuf'-(b*KfuiKuuKuu))*(iKuuKuf*b') + 2.*sum(sum(L'.*(L'*DKuf'*iKuuKuf))) - ...
                               sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
            
            for kk=1:length(ind)
              if savememory
                DKff=gpcf.fh.cfg(gpcf, x(ind{kk},:),[],[],i2);
              else
                DKff=DKffc{kk}{i2};
              end
              gdata(i1) = gdata(i1) ...
                  + 0.5.*(-b(ind{kk})*DKff*b(ind{kk})' ...
                          + 2.*b(ind{kk})*DKuf(:,ind{kk})'*iKuuKuf(:,ind{kk})*b(ind{kk})'- ...
                          b(ind{kk})*KfuiKuuKuu(ind{kk},:)*iKuuKuf(:,ind{kk})*b(ind{kk})' ...
                          + trace(La{kk}\DKff)...
                          - trace(L(ind{kk},:)*(L(ind{kk},:)'*DKff)) ...
                          + 2.*sum(sum(L(ind{kk},:)'.*(L(ind{kk},:)'*DKuf(:,ind{kk})'*iKuuKuf(:,ind{kk})))) - ...
                          sum(sum(L(ind{kk},:)'.*((L(ind{kk},:)'*KfuiKuuKuu(ind{kk},:))*iKuuKuf(:,ind{kk})))));                
            end
            gprior(i1) = gprior_cf(i2);
          end
          
          % Set the gradients of hyperparameter
          if length(gprior_cf) > np
            for i2=np+1:length(gprior_cf)
              i1 = i1+1;
              gdata(i1) = 0;
              gprior(i1) = gprior_cf(i2);
            end
          end
        end

      end
      
      % =================================================================
      % Gradient with respect to likelihood function parameters
      
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh, 'siteDeriv')

        [Ef, Varf] = gpep_pred(gp, x, y, x, 'tstind', gp.tr_index, 'z', z);
        
        sigm2_i = (Varf.^-1 - tautilde).^-1;
        mu_i = sigm2_i.*(Ef./Varf - nutilde);
        
        gdata_lik = 0;
        lik = gp.lik;
        for k1 = 1:length(y)
          gdata_lik = gdata_lik - lik.fh.siteDeriv(lik, y, k1, sigm2_i(k1), mu_i(k1), z);
        end

        % evaluate prior contribution for the gradient
        if isfield(gp.lik, 'p')
          g_logPrior = -lik.fh.lpg(lik);
        else
          g_logPrior = zeros(size(gdata_lik));
        end
        % set the gradients into vectors that will be returned
        gdata = [gdata gdata_lik];
        gprior = [gprior g_logPrior];
        i1 = length(gdata);
      end
      
      % =================================================================
      % Gradient with respect to inducing inputs
      
      if ~isempty(strfind(gp.infer_params, 'inducing'))
        if isfield(gp.p, 'X_u') && ~isempty(gp.p.X_u)
          m = size(gp.X_u,2);
          
          st=0;
          if ~isempty(gprior)
            st = length(gprior);
          end
          gdata(st+1:st+length(gp.X_u(:))) = 0;
          
          i1 = st+1;
          for i = 1:size(gp.X_u,1)
            if iscell(gp.p.X_u) % Own prior for each inducing input
              pr = gp.p.X_u{i};
              gprior(i1:i1+m) = pr.fh.lpg(gp.X_u(i,:), pr);
            else % One prior for all inducing inputs
              gprior(i1:i1+m-1) = gp.p.X_u.fh.lpg(gp.X_u(i,:), gp.p.X_u);
            end
            i1 = i1 + m;
          end
          
          % Loop over the  covariance functions
          for i=1:ncf            
            i1=st;
            gpcf = gp.cf{i};
            if savememory
              % If savememory option is used, just get the number of
              % covariates in X and calculate gradients later
              np=gpcf.fh.ginput(gpcf,u,[],0);
            else
              DKuu = gpcf.fh.ginput(gpcf, u);
              DKuf = gpcf.fh.ginput(gpcf, u, x);
              np=1;
            end
            
            for i3 = 1:np
              if savememory
                DKuu=gpcf.fh.ginput(gpcf,u,[],i3);
                DKuf=gpcf.fh.ginput(gpcf,u,x,i3);
              end
              for i2 = 1:length(DKuu)
                i1 = i1+1;
                KfuiKuuDKuu_u = iKuuKuf'*DKuu{i2};
              
                gdata(i1) = gdata(i1) - 0.5.*((2*b*DKuf{i2}'-(b*KfuiKuuDKuu_u))*(iKuuKuf*b') + 2.*sum(sum(L'.*((L'*DKuf{i2}')*iKuuKuf))) - ...
                                            sum(sum(L'.*((L'*KfuiKuuDKuu_u)*iKuuKuf))));
              
                for kk=1:length(ind)
                  gdata(i1) = gdata(i1) + 0.5.*(2.*b(ind{kk})*DKuf{i2}(:,ind{kk})'*iKuuKuf(:,ind{kk})*b(ind{kk})'- ...
                                              b(ind{kk})*KfuiKuuDKuu_u(ind{kk},:)*iKuuKuf(:,ind{kk})*b(ind{kk})' ...
                                              + 2.*sum(sum(L(ind{kk},:)'.*(L(ind{kk},:)'*DKuf{i2}(:,ind{kk})'*iKuuKuf(:,ind{kk})))) - ...
                                              sum(sum(L(ind{kk},:)'.*((L(ind{kk},:)'*KfuiKuuDKuu_u(ind{kk},:))*iKuuKuf(:,ind{kk})))));
                end
              end
            end
          end
        end
      end
      

    case {'CS+FIC'}
      % ============================================================
      % CS+FIC
      % ============================================================
      g_ind = zeros(1,numel(gp.X_u));
      gdata_ind = zeros(1,numel(gp.X_u));
      gprior_ind = zeros(1,numel(gp.X_u));

      u = gp.X_u;
      DKuu_u = 0;
      DKuf_u = 0;

      [e, edata, eprior, tautilde, nutilde, L, La, b] = gpep_e(w, gp, x, y, 'z', z);

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
      K_fu = gp_cov(gp, x, u);         % f x u
      K_uu = gp_trcov(gp, u);    % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu
      iKuuKuf = K_uu\K_fu';
      gp.cf = cf_orig;
      
      LD = ldlchol(La);
      siLa = spinv(LD,1);
      idiagLa = diag(siLa);
      LL = sum(L.*L,2);
      
      % =================================================================
      % Gradient with respect to covariance function parameters
      if ~isempty(strfind(gp.infer_params, 'covariance'))
        for i=1:ncf            
          i1=0;
          if ~isempty(gprior)
            i1 = length(gprior);
          end
          
          gpcf = gp.cf{i};
          
          % Evaluate the gradient for FIC covariance functions
          if ~isfield(gpcf,'cs')
            % Get the gradients of the covariance matrices 
            % and gprior from gpcf_* structures
            if savememory
              % If savememory option is used, just get the number of
              % hyperparameters and calculate gradients later
              np=gpcf.fh.cfg(gpcf,[],[],[],0);
            else
              DKffc = gpcf.fh.cfg(gpcf, x, [], 1);
              DKuuc = gpcf.fh.cfg(gpcf, u);
              DKufc = gpcf.fh.cfg(gpcf, u, x);
              np=length(DKuuc);
            end
            gprior_cf = -gpcf.fh.lpg(gpcf);

            for i2 = 1:np
              i1 = i1+1;
              if savememory
                DKff=gpcf.fh.cfg(gpcf,x,[],1,i2);
                DKuu=gpcf.fh.cfg(gpcf,u,[],[],i2);
                DKuf=gpcf.fh.cfg(gpcf,u,x,[],i2);
              else
                DKff=DKffc{i2};
                DKuu=DKuuc{i2};
                DKuf=DKufc{i2};
              end

              KfuiKuuKuu = iKuuKuf'*DKuu;
              gdata(i1) = -0.5.*((2*b*DKuf'-(b*KfuiKuuKuu))*(iKuuKuf*b') + 2.*sum(sum(L'.*(L'*DKuf'*iKuuKuf))) - ...
                                 sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
              
              gdata(i1) = gdata(i1) - 0.5.*(b.*DKff')*b';
              gdata(i1) = gdata(i1) + 0.5.*(2.*b.*sum(DKuf'.*iKuuKuf',2)'*b'- b.*sum(KfuiKuuKuu.*iKuuKuf',2)'*b');
              gdata(i1) = gdata(i1) + 0.5.*sum(idiagLa.*DKff - LL.*DKff);   % corrected
              gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf'.*iKuuKuf',2)) - sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
              
              %gdata(i1) = gdata(i1) + 0.5.*sum(sum(La\((2.*K_uf') - KfuiKuuKuu).*iKuuKuf',2));
              gdata(i1) = gdata(i1) + 0.5.*sum(sum(ldlsolve(LD,2.*DKuf' - KfuiKuuKuu).*iKuuKuf',2));
              gdata(i1) = gdata(i1) - 0.5.*( idiagLa'*(sum((2.*DKuf' - KfuiKuuKuu).*iKuuKuf',2)) ); % corrected                
              gprior(i1) = gprior_cf(i2);                    
            end                        

            % Evaluate the gradient for compact support covariance functions
          else
            % Get the gradients of the covariance matrices 
            % and gprior from gpcf_* structures
            if savememory
              % If savememory option is used, just get the number of
              % hyperparameters and calculate gradients later
              np=gpcf.fh.cfg(gpcf,[],[],[],0);
            else
              DKffc = gpcf.fh.cfg(gpcf, x);
              np=length(DKffc);
            end
            gprior_cf = -gpcf.fh.lpg(gpcf);
            for i2 = 1:np
              if savememory
                DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
              else
                DKff=DKffc{i2};
              end
              i1 = i1+1;
              gdata(i1) = 0.5*(sum(sum(siLa.*DKff',2)) - sum(sum(L.*(L'*DKff')')) - b*DKff*b');
              gprior(i1) = gprior_cf(i2);
            end
          end
          
          % Set the gradients of hyperparameter
          if length(gprior_cf) > np
            for i2=np+1:length(gprior_cf)
              i1 = i1+1;
              gdata(i1) = 0;
              gprior(i1) = gprior_cf(i2);
            end
          end
        end

      end         
      
      % =================================================================
      % Gradient with respect to likelihood function parameters
      
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh, 'siteDeriv')
        [Ef, Varf] = gpep_pred(gp, x, y, x, 'tstind', 1:n, 'z', z);
        sigm2_i = (Varf.^-1 - tautilde).^-1;
        mu_i = sigm2_i.*(Ef./Varf - nutilde);
        
        gdata_lik = 0;
        lik = gp.lik;
        for k1 = 1:length(y)
          gdata_lik = gdata_lik - lik.fh.siteDeriv(lik, y, k1, sigm2_i(k1), mu_i(k1), z);
        end
        % evaluate prior contribution for the gradient
        if isfield(gp.lik, 'p')
          g_logPrior = -lik.fh.lpg(lik);
        else
          g_logPrior = zeros(size(gdata_lik));
        end
        % set the gradients into vectors that will be returned
        gdata = [gdata gdata_lik];
        gprior = [gprior g_logPrior];
        i1 = length(gdata);
      end
      
      % =================================================================
      % Gradient with respect to inducing inputs
      
      if ~isempty(strfind(gp.infer_params, 'inducing'))
        if isfield(gp.p, 'X_u') && ~isempty(gp.p.X_u)
          m = size(gp.X_u,2);
          st=0;
          if ~isempty(gprior)
            st = length(gprior);
          end
          
          gdata(st+1:st+length(gp.X_u(:))) = 0;
          i1 = st+1;
          for i = 1:size(gp.X_u,1)
            if iscell(gp.p.X_u) % Own prior for each inducing input
              pr = gp.p.X_u{i};
              gprior(i1:i1+m) = pr.fh.lpg(gp.X_u(i,:), pr);
            else % One prior for all inducing inputs
              gprior(i1:i1+m-1) = gp.p.X_u.fh.lpg(gp.X_u(i,:), gp.p.X_u);
            end
            i1 = i1 + m;
          end
          
          for i=1:ncf
            i1=st;
            
            gpcf = gp.cf{i};            
            if ~isfield(gpcf,'cs')
              if savememory
                % If savememory option is used, just get the number of
                % covariates in X and calculate gradients later
                np=gpcf.fh.ginput(gpcf,u,[],0);
              else
                DKuu = gpcf.fh.ginput(gpcf, u);
                DKuf = gpcf.fh.ginput(gpcf, u, x);
                np=1;
              end
              
              for i3 = 1:np
                if savememory
                  DKuu=gpcf.fh.ginput(gpcf,u,[],i3);
                  DKuf=gpcf.fh.ginput(gpcf,u,x,i3);
                end
                for i2 = 1:length(DKuu)
                  i1 = i1 + 1;
                  KfuiKuuKuu = iKuuKuf'*DKuu{i2};
                
                  gdata(i1) = gdata(i1) - 0.5.*((2*b*DKuf{i2}'-(b*KfuiKuuKuu))*(iKuuKuf*b') + ...
                                              2.*sum(sum(L'.*(L'*DKuf{i2}'*iKuuKuf))) - sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
                  gdata(i1) = gdata(i1) + 0.5.*(2.*b.*sum(DKuf{i2}'.*iKuuKuf',2)'*b'- b.*sum(KfuiKuuKuu.*iKuuKuf',2)'*b');
                  gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf{i2}'.*iKuuKuf',2)) - ...
                                              sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
                
                  gdata(i1) = gdata(i1) + 0.5.*sum(sum(ldlsolve(LD,(2.*DKuf{i2}') - KfuiKuuKuu).*iKuuKuf',2));
                  gdata(i1) = gdata(i1) - 0.5.*( idiagLa'*(sum((2.*DKuf{i2}' - KfuiKuuKuu).*iKuuKuf',2)) ); % corrected
                  gprior(i1) = gprior_ind(i2);
                end
              end
            end
          end
        end
      end
      
    case {'DTC' 'SOR'}
      % ============================================================
      % DTC/SOR
      % ============================================================        
      g_ind = zeros(1,numel(gp.X_u));
      gdata_ind = zeros(1,numel(gp.X_u));
      gprior_ind = zeros(1,numel(gp.X_u));
      
      u = gp.X_u;
      DKuu_u = 0;
      DKuf_u = 0;

      [e, edata, eprior, tautilde, nutilde, L, La, b] = gpep_e(w, gp, x, y, 'z', z);

      K_fu = gp_cov(gp, x, u);         % f x u
      K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu        
      iKuuKuf = K_uu\K_fu';
      
      LL = sum(L.*L,2);
      iLav=1./La;
      
      % =================================================================
      % Gradient with respect to covariance function parameters
      if ~isempty(strfind(gp.infer_params, 'covariance'))
        for i=1:ncf            
          i1=0;
          if ~isempty(gprior)
            i1 = length(gprior);
          end
          
          gpcf = gp.cf{i};
          % Get the gradients of the covariance matrices 
          % and gprior from gpcf_* structures
          gpcf = gp.cf{i};
          if savememory
            % If savememory option is used, just get the number of
            % hyperparameters and calculate gradients later
            np=gpcf.fh.cfg(gpcf,[],[],[],0);
          else
%             DKffc = gpcf.fh.cfg(gpcf, x, [], 1);
            DKuuc = gpcf.fh.cfg(gpcf, u);
            DKufc = gpcf.fh.cfg(gpcf, u, x);
            np=length(DKuuc);
          end
          gprior_cf = -gpcf.fh.lpg(gpcf);
          
          for i2 = 1:np
            i1 = i1+1;
            if savememory
%               DKff=gpcf.fh.cfg(gpcf,x,[],1,i2);
              DKuu=gpcf.fh.cfg(gpcf,u,[],[],i2);
              DKuf=gpcf.fh.cfg(gpcf,u,x,[],i2);
            else
%               DKff=DKffc{i2};
              DKuu=DKuuc{i2};
              DKuf=DKufc{i2};
            end
            
            KfuiKuuKuu = iKuuKuf'*DKuu;
            gdata(i1) = -0.5.*((2*b*DKuf'-(b*KfuiKuuKuu))*(iKuuKuf*b'));
            gdata(i1) = gdata(i1) + 0.5.*(2.*(sum(iLav'*sum(DKuf'.*iKuuKuf',2))-sum(sum(L'.*(L'*DKuf'*iKuuKuf))))...
                                          - sum(iLav'*sum(KfuiKuuKuu.*iKuuKuf',2))+ sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
            gprior(i1) = gprior_cf(i2);
          end
          
          % Set the gradients of hyperparameter
          if length(gprior_cf) > np
            for i2=np+1:length(gprior_cf)
              i1 = i1+1;
              gdata(i1) = 0;
              gprior(i1) = gprior_cf(i2);
            end
          end
        end
        
      end
      
      % =================================================================
      % Gradient with respect to likelihood function parameters
      
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh, 'siteDeriv')
        [Ef, Varf] = gpep_pred(gp, x, y, x, 'tstind', 1:n, 'z', z);
        gdata_lik = 0;
        lik = gp.lik;
        for k1 = 1:length(y)
          sigm2_i = Varf(k1) ;
          mu_i = Ef(k1);
          gdata_lik = gdata_lik - lik.fh.siteDeriv(lik, y, k1, sigm2_i, mu_i, z);
        end
        % evaluate prior contribution for the gradient
        if isfield(gp.lik, 'p')
          g_logPrior = -lik.fh.lpg(lik);
        else
          g_logPrior = zeros(size(gdata_lik));
        end
        % set the gradients into vectors that will be returned
        gdata = [gdata gdata_lik];
        gprior = [gprior g_logPrior];
        i1 = length(gdata);
      end
      
      
      % =================================================================
      % Gradient with respect to inducing inputs
      
      if ~isempty(strfind(gp.infer_params, 'inducing'))
        if isfield(gp.p, 'X_u') && ~isempty(gp.p.X_u)
          m = size(gp.X_u,2);
          st=0;
          if ~isempty(gprior)
            st = length(gprior);
          end
          
          gdata(st+1:st+length(gp.X_u(:))) = 0;
          i1 = st+1;
          for i = 1:size(gp.X_u,1)
            if iscell(gp.p.X_u) % Own prior for each inducing input
              pr = gp.p.X_u{i};
              gprior(i1:i1+m) = pr.fh.lpg(gp.X_u(i,:), pr);
            else % One prior for all inducing inputs
              gprior(i1:i1+m-1) = gp.p.X_u.fh.lpg(gp.X_u(i,:), gp.p.X_u);
            end
            i1 = i1 + m;
          end
          
          for i=1:ncf
            i1=st;
            
            gpcf = gp.cf{i};
            if savememory
              % If savememory option is used, just get the number of
              % covariates in X and calculate gradients later
              np=gpcf.fh.ginput(gpcf,u,[],0); 
            else
              DKuu = gpcf.fh.ginput(gpcf, u);
              DKuf = gpcf.fh.ginput(gpcf, u, x);
              np=1;
            end
            
            for i3 = 1:np
              if savememory
                DKuu=gpcf.fh.ginput(gpcf,u,[],i3);
                DKuf=gpcf.fh.ginput(gpcf,u,x,i3);
              end
              for i2 = 1:length(DKuu)
                i1 = i1+1;
                
                KfuiKuuKuu = iKuuKuf'*DKuu{i2};
                gdata(i1) = gdata(i1) - 0.5.*((2*b*DKuf{i2}'-(b*KfuiKuuKuu))*(iKuuKuf*b'));
                gdata(i1) = gdata(i1) + 0.5.*(2.*(sum(iLav'*sum(DKuf{i2}'.*iKuuKuf',2))-sum(sum(L'.*(L'*DKuf{i2}'*iKuuKuf))))...
                                            - sum(iLav'*sum(KfuiKuuKuu.*iKuuKuf',2))+ sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));

                if strcmp(gp.type, 'VAR')
                  gdata(i1) = gdata(i1) + 0.5.*(0-2.*sum(iLav'*sum(DKuf{i2}'.*iKuuKuf',2)) + ...
                                              sum(iLav'*sum(KfuiKuuKuu.*iKuuKuf',2)));
                end
              end
            end
          end
        end
      end
      
      
    case 'VAR'
      % ============================================================
      % VAR
      % ============================================================        
      % NOTE! Not properly implemented as no analytical result has been
      % derived. Not suitable for large data sets.
      epsilon = 1.0e-6;
      
      func = fcnchk(@gpep_e, 3);
      %grad = fcnchk(grad, length(varargin));
      
      % Treat
      nparams = length(w);
      deltaf = zeros(1, nparams);
      step = zeros(1, nparams);
      for i = 1:nparams
        % Move a small way in the ith coordinate of w
        step(i) = 1.0;
        func = fcnchk(func, 3);
        fplus = func(w+epsilon.*step, gp,x,y);
        fminus = func(w-epsilon.*step, gp,x,y);
        %   fplus  = feval('linef_test', epsilon, func, w, step, varargin{:});
        %   fminus = feval('linef_test', -epsilon, func, w, step, varargin{:});
        % Use central difference formula for approximation
        deltaf(i) = 0.5*(fplus - fminus)/epsilon;
        step(i) = 0.0;
      end
      % gradient = feval(grad, w, varargin{:});
      % fprintf(1, 'Checking gradient ...\n\n');
      % fprintf(1, '   analytic   diffs     delta\n\n');
      % disp([gradient', deltaf', gradient' - deltaf'])
      
      %delta = gradient' - deltaf';
      gdata=deltaf;
      
      %gdata=numgrad_test(gp_pak(gp), @gpep_e, gp, x, y);
      gprior=0;

    case 'SSGP'        
      % ============================================================
      % SSGP
      % ============================================================        
      
      [e, edata, eprior, tautilde, nutilde, L, S, b] = gpep_e(w, gp, x, y, 'z', z);

      Phi = gp_trcov(gp, x);         % f x u
      m = size(Phi,2);

      SPhi = repmat(S,1,m).*Phi;                        
      % =================================================================
      % Evaluate the gradients from covariance functions
      for i=1:ncf            
        i1=0;
        if ~isempty(gprior)
          i1 = length(gprior);
        end
        
        gpcf = gp.cf{i};
        gpcf.GPtype = gp.type;
        % Covariance function parameters
        %--------------------------------------
        if ~isempty(strfind(gp.infer_params, 'covariance'))
          % Get the gradients of the covariance matrices 
          % and gprior from gpcf_* structures
          DKff = gpcf.fh.cfg(gpcf, x);
          gprior = -gpcf.fh.lpg(gpcf);
          i1 = i1+1;
          i2 = 1;
          
          % Evaluate the gradient with respect to magnSigma
          SDPhi = repmat(S,1,m).*DKff{i2};                
          
          gdata(i1) = 0.5*( sum(sum(SDPhi.*Phi,2)) + sum(sum(SPhi.*DKff{i2},2)) );
          gdata(i1) = gdata(i1) - 0.5*( sum(sum(L'.*(L'*DKff{i2}*Phi' + L'*Phi*DKff{i2}'),1)) );
          gdata(i1) = gdata(i1) - 0.5*(b*DKff{i2}*Phi' + b*Phi*DKff{i2}')*b';
          
          
          if isfield(gpcf.p.lengthScale, 'p') && ~isempty(gpcf.p.lengthScale.p)
            i1 = i1+1;
            if any(strcmp(fieldnames(gpcf.p.lengthScale.p),'nu'))
              i1 = i1+1;
            end
          end

          % Evaluate the gradient with respect to lengthScale
          for i2 = 2:length(DKff)
            i1 = i1+1;
            SDPhi = repmat(S,1,m).*DKff{i2};
            
            gdata(i1) = 0.5*( sum(sum(SDPhi.*Phi,2)) + sum(sum(SPhi.*DKff{i2},2)) );
            gdata(i1) = gdata(i1) - 0.5*( sum(sum(L'.*(L'*DKff{i2}*Phi' + L'*Phi*DKff{i2}'),1)) );
            gdata(i1) = gdata(i1) - 0.5*(b*DKff{i2}*Phi' + b*Phi*DKff{i2}')*b';
          end
        end
      end

      % likelihood parameters
      %--------------------------------------
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh, 'siteDeriv')
        [Ef, Varf] = gpep_pred(gp, x, y, x, param);
        gdata_lik = 0;
        lik = gp.lik;
        for k1 = 1:length(y)
          sigm2_i = Varf(k1) ;
          mu_i = Ef(k1);
          gdata_lik = gdata_lik + lik.fh.siteDeriv(lik, y, k1, sigm2_i, mu_i, z);
        end
        % evaluate prior contribution for the gradient
        if isfield(gp.lik, 'p')
          g_logPrior = -lik.fh.lpg(lik);
        else
          g_logPrior = zeros(size(gdata_lik));
        end
        % set the gradients into vectors that will be returned
        gdata = [gdata gdata_lik];
        gprior = [gprior g_logPrior];
        i1 = length(gdata);
      end
      
      
  end
  
  g = gdata + gprior;
  
end

