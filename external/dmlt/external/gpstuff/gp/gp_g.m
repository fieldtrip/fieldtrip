function [g, gdata, gprior] = gp_g(w, gp, x, y, varargin)
%GP_G  Evaluate the gradient of energy (GP_E) for Gaussian Process
%
%  Description
%    G = GP_G(W, GP, X, Y, OPTIONS) takes a full GP parameter
%    vector W, GP structure GP, a matrix X of input vectors and a
%    matrix Y of target vectors, and evaluates the gradient G of
%    the energy function (gp_e). Each row of X corresponds to one
%    input vector and each row of Y corresponds to one target
%    vector.
%
%    [G, GDATA, GPRIOR] = GP_G(W, GP, X, Y, OPTIONS) also returns
%    separately the data and prior contributions to the gradient.
%
%    OPTIONS is optional parameter-value pair
%      z - optional observed quantity in triplet (x_i,y_i,z_i)
%          Some likelihoods may use this. For example, in case of
%          Poisson likelihood we have z_i=E_i, that is, expected
%          value for ith case.
%
%  See also
%    GP_E, GP_PAK, GP_UNPAK, GPCF_*
%

% Copyright (c) 2007-2011 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari
% Copyright (c) 2010 Heikki Peura

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

if isfield(gp,'latent_method') && ~strcmp(gp.latent_method,'MCMC')
  % use an inference specific method
  fh_g = gp.fh.g;
  switch nargout 
    case {0 1}
      [g] = fh_g(w, gp, x, y, varargin{:});
    case 2
      [g, gdata] = fh_g(w, gp, x, y, varargin{:});
    case 3
      [g, gdata, gprior] = fh_g(w, gp, x, y, varargin{:});
  end
  return
end

ip=inputParser;
ip.FunctionName = 'GP_G';
ip.addRequired('w', @(x) isvector(x) && isreal(x));
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.parse(w, gp, x, y, varargin{:});
z=ip.Results.z;
if ~all(isfinite(w(:)));
  % instead of stopping to error, return NaN
  g=NaN;
  gdata = NaN;
  gprior = NaN;
  return;
end

% unpak the parameters
gp=gp_unpak(gp, w);
ncf = length(gp.cf);
if isfield(gp.lik, 'nondiagW')
  % Likelihoods with non-diagonal Hessian
  switch gp.lik.type
    case {'LGP', 'LGPC'}
      % Do nothing
    case {'Softmax', 'Multinom'}
      n=size(x,1);
      nout=size(y(:),1)./n;
%       [n,nout]=size(y);
      nl=cumsum([0 repmat(n,1,nout)]);
    otherwise
      n=size(x,1);
      nout=length(gp.comp_cf);
      
      % Help indices for latent processes
      if ~isfield(gp.lik, 'xtime')
        nl=[0 repmat(n,1,nout)];
      else
        xtime=gp.lik.xtime;
        ntime=size(xtime,1);
        nl=[0 ntime n];
      end
      nl=cumsum(nl);
  end
  if isfield(gp, 'comp_cf')  % own covariance for each ouput component
    multicf = true;
    if length(gp.comp_cf) ~= nout
      error('GP2_G: the number of component vectors in gp.comp_cf must be the same as number of outputs.')
    end
  else
    multicf = false;
  end
else
  n=size(x,1);
end

g = [];
gdata = [];
gprior = [];

if isfield(gp,'savememory') && gp.savememory
  savememory=1;
else
  savememory=0;
end

switch gp.type
  case 'FULL'
    % ============================================================
    % FULL
    % ============================================================
    
    if ~isfield(gp.lik, 'nondiagW') || ismember(gp.lik.type, {'LGP' 'LGPC'})
      % Evaluate covariance
      [K, C] = gp_trcov(gp,x);
      
      
      if issparse(C)
        % evaluate the sparse inverse
        [LD, notpositivedefinite] = ldlchol(C);
        invC = spinv(LD,1);
        if notpositivedefinite
          % instead of stopping to chol error, return NaN
          g=NaN;
          gdata = NaN;
          gprior = NaN;
          return;
        end
        if  ~isfield(gp,'meanf')
          b = ldlsolve(LD,y);
        else
          [invNM invAt HinvC]=mean_gf(gp,x,C,LD,[],[],y,'gaussian');
        end
      else
        % evaluate the full inverse
        ws1=warning('off','MATLAB:nearlySingularMatrix');
        ws2=warning('off','MATLAB:SingularMatrix');
        invC = inv(C);
        if  ~isfield(gp,'meanf')
          b = C\y;
        else
          [invNM invAt HinvC]=mean_gf(gp,x,C,invC,[],[],y,'gaussian');
        end
        warning(ws1);
        warning(ws2);
      end
    else
      b = zeros(nl(end),1);
      y=y(:);
      
      switch gp.lik.type
        case 'Coxph'
          % In Cox-Ph, latent processes have different inputs so stacking in invC
          % is not possible.
          invC = zeros(nl(end),nl(end));
          if multicf
            [tmp, C] = gp_trcov(gp, xtime, gp.comp_cf{1});
            invC(1+nl(1):nl(2),1+nl(1):nl(2)) = inv(C);
            b(nl(1)+1:nl(2)) = C\y(nl(1)+1:nl(2));
            [tmp, C] = gp_trcov(gp, x, gp.comp_cf{2});
            invC(1+nl(2):nl(3),1+nl(2):nl(3)) = inv(C);
            b(1+nl(2):nl(3)) = C\y(1+nl(2):nl(3));
          else
            error('Specify covariance function for time process and input process, when using Cox-Ph likelihood');
          end
        otherwise
          invC = zeros(n,n,nout);
          if multicf
            for i1=1:nout
              [tmp, C] = gp_trcov(gp, x, gp.comp_cf{i1});
              invC(:,:,i1) = inv(C);
              b(1+nl(i1):nl(i1+1)) = C\y(1+nl(i1):nl(i1+1));
              %             b(:,i1) = C\y(:,i1);
            end
          else
            [tmp, C] = gp_trcov(gp, x);
            invCtmp = inv(C);
            for i1=1:nout
              invC(:,:,i1) = invCtmp;
              b(1+nl(i1):nl(i1+1)) = C\y(1+nl(i1):nl(i1+1));
              %             b(:,i1) = C\y(:,i1);
            end
          end
      end
    end

    % =================================================================
    % Gradient with respect to covariance function parameters
    i1=0;
    if ~isempty(strfind(gp.infer_params, 'covariance'))
      for i=1:ncf
        if ~isempty(gprior)
          i1 = length(gprior);
        end
        
        gpcf = gp.cf{i};
        
        if isfield(gp.lik, 'nondiagW') && ~ismember(gp.lik.type, {'LGP' 'LGPC'})
          % check in which components the covariance function is present
          % for likelihoods with non-diagonal Hessian
          do = false(nout,1);
          if multicf
            for z1=1:nout
              if any(gp.comp_cf{z1}==i)
                do(z1) = true;
              end
            end
          else
            do = true(nout,1);
          end         
        end
        
        if ~(isfield(gp,'derivobs') && gp.derivobs)
          % No derivative observations
          if ~savememory
            if isfield(gp.lik, 'nondiagW') && isfield(gp,'comp_cf') && ~isempty(intersect(gp.comp_cf{1},i)) && isfield(gp.lik, 'xtime')
              DKffc = gpcf.fh.cfg(gpcf, xtime);
            else
              DKffc = gpcf.fh.cfg(gpcf, x);
            end
            np=length(DKffc);
          else
            % If savememory option is used, just get the number of
            % hyperparameters and calculate gradients later
            np=gpcf.fh.cfg(gpcf,[],[],[],0);
          end
          gprior_cf = -gpcf.fh.lpg(gpcf);
        else
          [n m]=size(x);
          %Case: input dimension is 1
          if m==1

            DKffa = gpcf.fh.cfg(gpcf, x);
            DKdf = gpcf.fh.cfdg(gpcf, x);
            DKdd = gpcf.fh.cfdg2(gpcf, x);
            gprior_cf = -gpcf.fh.lpg(gpcf);

            % DKff{1} -- d K / d magnSigma2
            % DKff{2} -- d K / d lengthScale
            DKffc{1} = [DKffa{1}, DKdf{1}'; DKdf{1}, DKdd{1}];
            DKffc{2} = [DKffa{2}, DKdf{2}'; DKdf{2}, DKdd{2}];
            np=2;
            
            %Case: input dimension is >1    
          else
            DKffa = gpcf.fh.cfg(gpcf, x);
            DKdf = gpcf.fh.cfdg(gpcf, x);
            DKdd = gpcf.fh.cfdg2(gpcf, x);
            gprior_cf = -gpcf.fh.lpg(gpcf);

            %Check whether ARD method is in use (with gpcf_sexp)
            Ard=length(gpcf.lengthScale);
            
            % DKff{1} - d K / d magnSigma2
            % DKff{2:end} - d K / d lengthScale(1:end)
            for i=1:2
              DKffc{i}=[DKffa{i} DKdf{i}';DKdf{i} DKdd{i}];
            end
            
            %If ARD is in use
            if Ard>1
              for i=2+1:2+Ard-1
                DKffc{i}=[DKffa{i} DKdf{i}';DKdf{i} DKdd{i}];
              end  
            end
            np=length(DKffc);
          end
        end
        
        % Are there specified mean functions
        if  ~isfield(gp,'meanf')
          % Evaluate the gradient with respect to covariance function
          % parameters
          for i2 = 1:np
            if savememory
              if isfield(gp.lik, 'nondiagW') && isfield(gp,'comp_cf') && ~isempty(intersect(gp.comp_cf{1},i)) && isfield(gp.lik, 'xtime')
                DKff=gpcf.fh.cfg(gpcf,xtime,[],[],i2);
              else
                DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
              end
            else
              DKff=DKffc{i2};
            end
            i1 = i1+1;
            if ~isfield(gp.lik, 'nondiagW')
              Bdl = b'*(DKff*b);
              Cdl = sum(sum(invC.*DKff)); % help arguments
            else
              % Non-diagonalizable likelihoods
              Bdl=0; Cdl=0;
              if isfield(gp.lik,'xtime');
                if do(1)
                  Bdl = Bdl + b(1:ntime)'*(DKff*b(1:ntime));
                  Cdl = Cdl + sum(sum(invC(1:ntime,1:ntime).*DKff)); % help arguments
                end
                if do(2)
                  Bdl = Bdl + b(ntime+1:end)'*(DKff*b(ntime+1:end));
                  Cdl = Cdl + sum(sum(invC(ntime+1:end,ntime+1:end).*DKff)); % help arguments
                end
              else
                for z1=1:nout
                  if do(z1)
                    Bdl = Bdl + b(1+nl(z1):nl(z1+1))'*(DKff*b(1+nl(z1):nl(z1+1)));
                    Cdl = Cdl + sum(sum(invC(:,:,z1).*DKff)); % help arguments
                  end
                end
              end
            end
            gdata(i1)=0.5.*(Cdl - Bdl);
            gprior(i1) = gprior_cf(i2);
          end
        else
          for i2 = 1:np
            if savememory
              DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
            else
              DKff=DKffc{i2};
            end
            i1=i1+1;
            dA = -1*HinvC*DKff*HinvC';                  % d A / d th
            trA = sum(invAt(:).*dA(:));                 % d log(|A|) / dth
            dMNM = invNM'*(DKff*invNM);           % d M'*N*M / d th
            trK = sum(sum(invC.*DKff));       % d log(Ky�?�) / d th
            gdata(i1)=0.5*(-1*dMNM + trK + trA);
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
    % Gradient with respect to Gaussian likelihood function parameters
    if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh,'trcov')
      % Evaluate the gradient from Gaussian likelihood
      DCff = gp.lik.fh.cfg(gp.lik, x);
      gprior_lik = -gp.lik.fh.lpg(gp.lik);
      for i2 = 1:length(DCff)
        i1 = i1+1;
        if ~isfield(gp,'meanf')
          if size(DCff{i2}) > 1
            yKy = b'*(DCff{i2}*b);
            trK = sum(sum(invC.*DCff{i2})); % help arguments
            gdata_zeromean(i1)=0.5.*(trK - yKy);
          else 
            yKy=DCff{i2}.*(b'*b);
            trK = DCff{i2}.*(trace(invC));
            gdata_zeromean(i1)=0.5.*(trK - yKy);
          end
          gdata(i1)=gdata_zeromean(i1);
        else
          if size(DCff{i2}) > 1
            trK = sum(sum(invC.*DCff{i2})); % help arguments
          else 
            trK = DCff{i2}.*(trace(invC));
          end
          dA = -1*HinvC*DCff{i2}*HinvC';            % d A / d th
          trA = sum(invAt(:).*dA(:));               % d log(|A|) / dth
          dMNM = invNM'*(DCff{i2}*invNM);           % d M'*N*M / d th
          gdata(i1)=0.5*(-1*dMNM + trA + trK);
        end
        gprior(i1) = gprior_lik(i2);
      end
      
      % Set the gradients of hyperparameter
      if length(gprior_lik) > length(DCff)
        for i2=length(DCff)+1:length(gprior_lik)
          i1 = i1+1;
          gdata(i1) = 0;
          gprior(i1) = gprior_lik(i2);
        end
      end
    end
    
    
    if ~isempty(strfind(gp.infer_params, 'mean')) && isfield(gp,'meanf')
        notpositivedefinite2 = 0; notpositivedefinite3 = 0;
        
        nmf=numel(gp.meanf);
        [H,b,B]=mean_prep(gp,x,[]);
        M = H'*b-y;
        
        if issparse(C)
            [LD, notpositivedefinite] = ldlchol(C);
            if ~notpositivedefinite
                KH = ldlsolve(LD, H');
            end
            [LB, notpositivedefinite2] = chol(B);
            if ~notpositivedefinite2
                A = LB\(LB'\eye(size(B))) + H*KH;
                LA = chol(A);
                a = ldlsolve(LD, M) - KH*(LA\(LA'\(KH'*M)));
                iNH = ldlsolve(LD, H') - KH*(LA\(LA'\(KH'*H')));
            end
        else
            N = C + H'*B*H;
            [LN, notpositivedefinite3] = chol(N);
            if ~notpositivedefinite3
                a = LN\(LN'\M);
                iNH = LN\(LN'\H');
            end
        end
        if (~notpositivedefinite2 && ~notpositivedefinite3)
            Ha=H*a;
            g_bb = (-H*a)';     % b and B parameters are log transformed in packing 
            indB = find(B>0);
            for i=1:length(indB)
                Bt = zeros(size(B)); Bt(indB(i))=1;
                BH = Bt*H;
                g_B(i) = 0.5* ( Ha'*Bt*Ha - sum(sum(iNH.*(BH'))) );
            end
            g_BB = g_B.*B(indB)';
            for i=1:nmf
                gpmf = gp.meanf{i};
                [lpg_b, lpg_B] = gpmf.fh.lpg(gpmf);
                ll=length(lpg_b);
                gdata = [gdata -g_bb((i-1)*ll+1:i*ll)];
                gprior = [gprior -lpg_b];
                ll=length(lpg_B);
                gdata = [gdata -g_B((i-1)*ll+1:i*ll)];
                gprior = [gprior -lpg_B];
            end
        else
          g=NaN;
          gdata = NaN;
          gprior = NaN;
          return
        end
    end
    
    g = gdata + gprior;
    
  case 'FIC'
    % ============================================================
    % FIC
    % ============================================================
    g_ind = zeros(1,numel(gp.X_u));
    gdata_ind = zeros(1,numel(gp.X_u));
    gprior_ind = zeros(1,numel(gp.X_u));

    u = gp.X_u;
    DKuu_u = 0;
    DKuf_u = 0;

    % First evaluate the needed covariance matrices
    % v defines that parameter is a vector
    [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
    K_fu = gp_cov(gp, x, u);         % f x u
    K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu
    [Luu, notpositivedefinite] = chol(K_uu,'lower');
    if notpositivedefinite
      % If not positive definite, return NaN
      g=NaN; gdata=NaN; gprior=NaN;
      return;
    end
    % Evaluate the Lambda (La)
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    % Here we need only the diag(Q_ff), which is evaluated below
    B=Luu\(K_fu');
    Qv_ff=sum(B.^2)';
    Lav = Cv_ff-Qv_ff;   % 1 x f, Vector of diagonal elements
                         % iLaKfu = diag(inv(Lav))*K_fu = inv(La)*K_fu
    iLaKfu = zeros(size(K_fu));  % f x u,
    for i=1:n
      iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
    end
    % ... then evaluate some help matrices.
    % A = K_uu+K_uf*inv(La)*K_fu
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;               % Ensure symmetry
    [A, notpositivedefinite] = chol(A,'upper');
    if notpositivedefinite
      % If not positive definite, return NaN
      g=NaN; gdata=NaN; gprior=NaN;
      return;
    end
    L = iLaKfu/A;
    b = y'./Lav' - (y'*L)*L';
    iKuuKuf = Luu'\(Luu\K_fu');
    La = Lav;
    LL = sum(L.*L,2);
    
    % =================================================================
    % Gradient with respect to covariance function parameters
    if ~isempty(strfind(gp.infer_params, 'covariance'))
      % Loop over the covariance functions
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
          DKffc = gpcf.fh.cfg(gpcf, x, [], 1);
          DKuuc = gpcf.fh.cfg(gpcf, u);
          DKufc = gpcf.fh.cfg(gpcf, u, x);
          np=length(DKffc);
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
          gdata(i1) = gdata(i1) + 0.5.*(sum(DKff./La) - sum(LL.*DKff));
          gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf'.*iKuuKuf',2)) - sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
          gprior(i1) = gprior_cf(i2);
        end
        
        % Set the gradients of hyperparameter
        if length(gprior_cf) > np
          for i2=length(DKff)+1:length(gprior_cf)
            i1 = i1+1;
            gdata(i1) = 0;
            gprior(i1) = gprior_cf(i2);
          end
        end
      end
    end

    % =================================================================
    % Gradient with respect to Gaussian likelihood function parameters
    if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh,'trcov')
      % Evaluate the gradient from Gaussian likelihood
      DCff = gp.lik.fh.cfg(gp.lik, x);
      gprior_lik = -gp.lik.fh.lpg(gp.lik);
      for i2 = 1:length(DCff)
        i1 = i1+1;
        gdata(i1)= -0.5*DCff{i2}.*b*b';
        gdata(i1)= gdata(i1) + 0.5*sum(DCff{i2}./La-sum(L.*L,2).*DCff{i2});
        gprior(i1) = gprior_lik(i2);
      end
      % Set the gradients of hyperparameter
      if length(gprior_lik) > length(DCff)
        for i2=length(DCff)+1:length(gprior_lik)
          i1 = i1+1;
          gdata(i1) = 0;
          gprior(i1) = gprior_lik(i2);
        end
      end               
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
            gprior(i1:i1+m) = -pr.fh.lpg(gp.X_u(i,:), pr);
          else % One prior for all inducing inputs
            gprior(i1:i1+m-1) = -gp.p.X_u.fh.lpg(gp.X_u(i,:), gp.p.X_u);
          end
          i1 = i1 + m;
        end
        
        % Loop over the covariance functions
        for i=1:ncf
          i1 = st;
          gpcf = gp.cf{i};
          
          if savememory
            % If savememory option is used, just get the number of
            % covariates in X and calculate gradients later
            np=gpcf.fh.ginput(gpcf,u,[],0);
          else
            np=1;
            DKuu = gpcf.fh.ginput(gpcf, u);
            DKuf = gpcf.fh.ginput(gpcf, u, x);
          end
          for i3=1:np
            if savememory
              DKuu=gpcf.fh.ginput(gpcf,u,[],i3);
              DKuf=gpcf.fh.ginput(gpcf,u,x,i3);
            end
            for i2 = 1:length(DKuu)
              i1=i1+1;
              KfuiKuuKuu = iKuuKuf'*DKuu{i2};
            
              gdata(i1) = gdata(i1) - 0.5.*((2*b*DKuf{i2}'-(b*KfuiKuuKuu))*(iKuuKuf*b') + ...
                                          2.*sum(sum(L'.*(L'*DKuf{i2}'*iKuuKuf))) - sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
              gdata(i1) = gdata(i1) + 0.5.*(2.*b.*sum(DKuf{i2}'.*iKuuKuf',2)'*b'- b.*sum(KfuiKuuKuu.*iKuuKuf',2)'*b');
              gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf{i2}'.*iKuuKuf',2)) - ...
                                          sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
            end
          end
        end
      end
    end
    
    g = gdata + gprior;
    
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

    % First evaluate the needed covariance matrices
    % if they are not in the memory
    % v defines that parameter is a vector
    K_fu = gp_cov(gp, x, u);         % f x u
    K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu
    [Luu, notpositivedefinite] = chol(K_uu,'lower');
    if notpositivedefinite
      % If not positive definite, return NaN
      g=NaN; gdata=NaN; gprior=NaN;
      return;
    end
    % Evaluate the Lambda (La)
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    % Here we need only the diag(Q_ff), which is evaluated below
    %B=K_fu/Luu;
    B=Luu\K_fu';
    iLaKfu = zeros(size(K_fu));  % f x u
    for i=1:length(ind)
      Qbl_ff = B(:,ind{i})'*B(:,ind{i});
      [Kbl_ff, Cbl_ff] = gp_trcov(gp, x(ind{i},:));
      la = Cbl_ff - Qbl_ff;
      La{i} = (la + la')./2;
      iLaKfu(ind{i},:) = La{i}\K_fu(ind{i},:);
    end
    % ... then evaluate some help matrices.
    % A = chol(K_uu+K_uf*inv(La)*K_fu))
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;            % Ensure symmetry

    [LA,notpositivedefinite]=chol(A,'upper');
    if notpositivedefinite
      % instead of stopping to chol error, return NaN
      g=NaN; gdata = NaN; gprior = NaN;
      return;
    end
    L = iLaKfu/LA;
    b = zeros(1,n);
    b_apu=(y'*L)*L';
    for i=1:length(ind)
      b(ind{i}) = y(ind{i})'/La{i} - b_apu(ind{i});
    end
    iKuuKuf = Luu'\(Luu\K_fu');
    
    % =================================================================
    % Gradient with respect to covariance function parameters

    if ~isempty(strfind(gp.infer_params, 'covariance'))
      % Loop over the  covariance functions
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
              DKff=gpcf.fh.cfg(gpcf,x(ind{kk},:),[],[],i2);
            else
              DKff=DKffc{kk}{i2};
            end
            gdata(i1) = gdata(i1) ...
                + 0.5.*(-b(ind{kk})*DKff*b(ind{kk})' ...
                        + 2.*b(ind{kk})*DKuf(:,ind{kk})'*iKuuKuf(:,ind{kk})*b(ind{kk})'- ...
                        b(ind{kk})*KfuiKuuKuu(ind{kk},:)*iKuuKuf(:,ind{kk})*b(ind{kk})' ... 
                        +trace(La{kk}\DKff)...                                
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
    % Gradient with respect to Gaussian likelihood function parameters
    if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh,'trcov')
      % Evaluate the gradient from Gaussian likelihood
      DCff = gp.lik.fh.cfg(gp.lik, x);
      gprior_lik = -gp.lik.fh.lpg(gp.lik);
      for i2 = 1:length(DCff)
        i1 = i1+1;
        gdata(i1)= -0.5*DCff{i2}.*b*b';            
        for kk=1:length(ind)
          gdata(i1)= gdata(i1) + 0.5*trace((inv(La{kk})-L(ind{kk},:)*L(ind{kk},:)')).*DCff{i2};
        end
        gprior(i1) = gprior_lik(i2);
      end
      % Set the gradients of hyperparameter
      if length(gprior_lik) > length(DCff)
        for i2=length(DCff)+1:length(gprior_lik)
          i1 = i1+1;
          gdata(i1) = 0;
          gprior(i1) = gprior_lik(i2);
        end
      end
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
            gprior(i1:i1+m) = -pr.fh.lpg(gp.X_u(i,:), pr);
          else % One prior for all inducing inputs
            gprior(i1:i1+m-1) = -gp.p.X_u.fh.lpg(gp.X_u(i,:), gp.p.X_u);
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
            np=1;
            DKuu = gpcf.fh.ginput(gpcf, u);
            DKuf = gpcf.fh.ginput(gpcf, u, x);
          end
          
          for i3 = 1:np
            if savememory
              DKuu=gpcf.fh.ginput(gpcf, u, [], i3);
              DKuf=gpcf.fh.ginput(gpcf, u, x, i3);
            end
            for i2 = 1:length(DKuu)
              i1 = i1+1;
              KfuiKuuDKuu_u = iKuuKuf'*DKuu{i2};
              gdata(i1) = gdata(i1) -0.5.*((2*b*DKuf{i2}'-(b*KfuiKuuDKuu_u))*(iKuuKuf*b') + 2.*sum(sum(L'.*((L'*DKuf{i2}')*iKuuKuf))) - ...
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
    g = gdata + gprior;
    
  case 'CS+FIC'
    % ============================================================
    % CS+FIC
    % ============================================================
    g_ind = zeros(1,numel(gp.X_u));
    gdata_ind = zeros(1,numel(gp.X_u));
    gprior_ind = zeros(1,numel(gp.X_u));

    u = gp.X_u;
    DKuu_u = 0;
    DKuf_u = 0;

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

    % First evaluate the needed covariance matrices
    % if they are not in the memory
    % v defines that parameter is a vector
    [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
    K_fu = gp_cov(gp, x, u);         % f x u
    K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu
    [Luu, notpositivedefinite] = chol(K_uu,'lower');
    if notpositivedefinite
      % If not positive definite, return NaN
      g=NaN; gdata=NaN; gprior=NaN;
      return;
    end
    % Evaluate the Lambda (La)
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    % Here we need only the diag(Q_ff), which is evaluated below
    B=Luu\(K_fu');
    Qv_ff=sum(B.^2)';
    Lav = Cv_ff-Qv_ff;   % 1 x f, Vector of diagonal elements

    gp.cf = cf2;
    K_cs = gp_trcov(gp,x);
    La = sparse(1:n,1:n,Lav,n,n) + K_cs;
    gp.cf = cf_orig;

    LD = ldlchol(La);
    %        iLaKfu = La\K_fu;
    iLaKfu = ldlsolve(LD, K_fu);

    % ... then evaluate some help matrices.
    % A = chol(K_uu+K_uf*inv(La)*K_fu))
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;            % Ensure symmetry
    [LA, notpositivedefinite] = chol(A,'upper');
    if notpositivedefinite
      % If not positive definite, return NaN
      g=NaN; gdata=NaN; gprior=NaN;
      return;
    end
    L = iLaKfu/LA;
    %b = y'/La - (y'*L)*L';
    b = ldlsolve(LD,y)' - (y'*L)*L';
    
    siLa = spinv(La);
    idiagLa = diag(siLa);
    iKuuKuf = K_uu\K_fu';
    LL = sum(L.*L,2);
    
    % =================================================================
    % Gradient with respect to covariance function parameters
    if ~isempty(strfind(gp.infer_params, 'covariance'))
      % Loop over covariance functions 
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
            
            temp1 = sum(KfuiKuuKuu.*iKuuKuf',2);
            temp2 = sum(DKuf'.*iKuuKuf',2);
            temp3 = 2.*DKuf' - KfuiKuuKuu;
            gdata(i1) = gdata(i1) - 0.5.*(b.*DKff')*b';
            gdata(i1) = gdata(i1) + 0.5.*(2.*b.*temp2'*b'- b.*temp1'*b');
            gdata(i1) = gdata(i1) + 0.5.*(sum(idiagLa.*DKff - LL.*DKff));   % corrected
            gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*temp2) - sum(LL.*temp1));
            
            %gdata(i1) = gdata(i1) + 0.5.*sum(sum(La\((2.*K_uf') - KfuiKuuKuu).*iKuuKuf',2));
            gdata(i1) = gdata(i1) + 0.5.*sum(sum(ldlsolve(LD,temp3).*iKuuKuf',2));
            gdata(i1) = gdata(i1) - 0.5.*( idiagLa'*(sum(temp3.*iKuuKuf',2)) ); % corrected                
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
            i1 = i1+1;
            if savememory
              DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
            else
              DKff=DKffc{i2};
            end
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
    % Gradient with respect to Gaussian likelihood function parameters
    if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh,'trcov')
      % Evaluate the gradient from Gaussian likelihood
      DCff = gp.lik.fh.cfg(gp.lik, x);
      gprior_lik = -gp.lik.fh.lpg(gp.lik);
      for i2 = 1:length(DCff)
        i1 = i1+1;
        gdata(i1)= -0.5*DCff{i2}.*b*b';
        gdata(i1)= gdata(i1) + 0.5*sum(idiagLa-LL).*DCff{i2};
        gprior(i1) = gprior_lik(i2);
      end
      
      % Set the gradients of hyperparameter
      if length(gprior_lik) > length(DCff)
        for i2=length(DCff)+1:length(gprior_lik)
          i1 = i1+1;
          gdata(i1) = 0;
          gprior(i1) = gprior_lik(i2);
        end
      end
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
            gprior(i1:i1+m) = -pr.fh.lpg(gp.X_u(i,:), pr);
          else % One prior for all inducing inputs
            gprior(i1:i1+m-1) = -gp.p.X_u.fh.lpg(gp.X_u(i,:), gp.p.X_u);
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
              np=1;
              DKuu = gpcf.fh.ginput(gpcf, u);
              DKuf = gpcf.fh.ginput(gpcf, u, x);
            end
            
            for i3 = 1:np
              if savememory
                DKuu = gpcf.fh.ginput(gpcf,u,[],i3);
                DKuf = gpcf.fh.ginput(gpcf,u,x,i3);
              end
              for i2 = 1:length(DKuu)
                i1 = i1+1;
                KfuiKuuKuu = iKuuKuf'*DKuu{i2};
              
                gdata(i1) = gdata(i1) - 0.5.*((2*b*DKuf{i2}'-(b*KfuiKuuKuu))*(iKuuKuf*b') + ...
                                            2.*sum(sum(L'.*(L'*DKuf{i2}'*iKuuKuf))) - sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
                gdata(i1) = gdata(i1) + 0.5.*(2.*b.*sum(DKuf{i2}'.*iKuuKuf',2)'*b'- b.*sum(KfuiKuuKuu.*iKuuKuf',2)'*b');
                gdata(i1) = gdata(i1) + 0.5.*(2.*sum(sum(L.*L,2).*sum(DKuf{i2}'.*iKuuKuf',2)) - ...
                                            sum(sum(L.*L,2).*sum(KfuiKuuKuu.*iKuuKuf',2)));
              
                gdata(i1) = gdata(i1) + 0.5.*sum(sum(ldlsolve(LD,(2.*DKuf{i2}') - KfuiKuuKuu).*iKuuKuf',2));
                gdata(i1) = gdata(i1) - 0.5.*( idiagLa'*(sum((2.*DKuf{i2}' - KfuiKuuKuu).*iKuuKuf',2)) ); % corrected
              end
            end
          end
        end
      end
    end
    
    g = gdata + gprior;
    
  case {'DTC' 'VAR' 'SOR'}
    % ============================================================
    % DTC/VAR/SOR
    % ============================================================
    g_ind = zeros(1,numel(gp.X_u));
    gdata_ind = zeros(1,numel(gp.X_u));
    gprior_ind = zeros(1,numel(gp.X_u));

    u = gp.X_u;
    DKuu_u = 0;
    DKuf_u = 0;

    % First evaluate the needed covariance matrices
    % v defines that parameter is a vector
    [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
    K_fu = gp_cov(gp, x, u);         % f x u
    K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu
    [Luu, notpositivedefinite] = chol(K_uu,'lower');
    if notpositivedefinite
      % If not positive definite, return NaN
      g=NaN; gdata=NaN; gprior=NaN;
      return;
    end
    % Evaluate the Lambda (La)
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    % Here we need only the diag(Q_ff), which is evaluated below
    B=Luu\(K_fu');
    Qv_ff=sum(B.^2)';
    Lav = Cv_ff-Kv_ff;   % 1 x f, Vector of diagonal elements
                         % iLaKfu = diag(inv(Lav))*K_fu = inv(La)*K_fu
    iLaKfu = zeros(size(K_fu));  % f x u,
    for i=1:n
      iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
    end
    % ... then evaluate some help matrices.
    % A = K_uu+K_uf*inv(La)*K_fu
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;               % Ensure symmetry
    [LA, notpositivedefinite] = chol(A);
    if notpositivedefinite
      % If not positive definite, return NaN
      g=NaN; gdata=NaN; gprior=NaN;
      return;
    end
    L = iLaKfu/LA;
    b = y'./Lav' - (y'*L)*L';
    iKuuKuf = Luu'\(Luu\K_fu');
    La = Lav;
    LL = sum(L.*L,2);
    iLav=1./Lav;
    
    LL1=iLav-LL;
    
    % =================================================================
    
    if ~isempty(strfind(gp.infer_params, 'covariance'))
      % Loop over the covariance functions
      for i=1:ncf            
        i1=0;
        if ~isempty(gprior)
          i1 = length(gprior);
        end
        
        % Get the gradients of the covariance matrices 
        % and gprior from gpcf_* structures
        gpcf = gp.cf{i};
        if savememory
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
            % If savememory option is used, just get the number of
            % hyperparameters and calculate gradients later
            DKff=gpcf.fh.cfg(gpcf,x,[],1,i2);
            DKuu=gpcf.fh.cfg(gpcf,u,[],[],i2);
            DKuf=gpcf.fh.cfg(gpcf,u,x,[],i2);
          else
            DKff=DKffc{i2};
            DKuu=DKuuc{i2};
            DKuf=DKufc{i2};
          end
          
          KfuiKuuKuu = iKuuKuf'*DKuu;
          gdata(i1) = -0.5.*((2*b*DKuf'-(b*KfuiKuuKuu))*(iKuuKuf*b'));
          gdata(i1) = gdata(i1) + 0.5.*(2.*(sum(iLav'*sum(DKuf'.*iKuuKuf',2))-sum(sum(L'.*(L'*DKuf'*iKuuKuf))))...
                                        - sum(iLav'*sum(KfuiKuuKuu.*iKuuKuf',2))+ sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
          
          if strcmp(gp.type, 'VAR')
            gdata(i1) = gdata(i1) + 0.5.*(sum(iLav.*DKff)-2.*sum(iLav'*sum(DKuf'.*iKuuKuf',2)) + ...
                                          sum(iLav'*sum(KfuiKuuKuu.*iKuuKuf',2))); % trace-term derivative
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
    % Gradient with respect to Gaussian likelihood function parameters
    if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh,'trcov')
      % Evaluate the gradient from Gaussian likelihood
      DCff = gp.lik.fh.cfg(gp.lik, x);
      gprior_lik = -gp.lik.fh.lpg(gp.lik);
      for i2 = 1:length(DCff)
        i1 = i1+1;
        gdata(i1)= -0.5*DCff{i2}.*b*b';
        gdata(i1)= gdata(i1) + 0.5*sum(DCff{i2}./La-sum(L.*L,2).*DCff{i2});
        if strcmp(gp.type, 'VAR')
          gdata(i1)= gdata(i1) + 0.5*(sum((Kv_ff-Qv_ff)./La));
        end
        
        gprior(i1) = gprior_lik(i2);                        
      end
      % Set the gradients of hyperparameter
      if length(gprior_lik) > length(DCff)
        for i2=length(DCff)+1:length(gprior_lik)
          i1 = i1+1;
          gdata(i1) = 0;
          gprior(i1) = gprior_lik(i2);
        end
      end               
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
            gprior(i1:i1+m) = -pr.fh.lpg(gp.X_u(i,:), pr);
          else % One prior for all inducing inputs
            gprior(i1:i1+m-1) = -gp.p.X_u.fh.lpg(gp.X_u(i,:), gp.p.X_u);
          end
          i1 = i1 + m;
        end
        
        % Loop over the covariance functions
        for i=1:ncf
          i1 = st;
          gpcf = gp.cf{i};
          if savememory
            % If savememory option is used, just get the number of
            % covariates in X and calculate gradients later
            np=gpcf.fh.ginput(gpcf,u,[],0);
          else
            np=1;
            DKuu = gpcf.fh.ginput(gpcf, u);
            DKuf = gpcf.fh.ginput(gpcf, u, x);
          end
          
          for i3 = 1:np
            if savememory
              DKuu=gpcf.fh.ginput(gpcf,u,[],i3);
              DKuf=gpcf.fh.ginput(gpcf,u,x,i3);
            end
            for i2 = 1:length(DKuu)
              i1=i1+1;
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
    
    g = gdata + gprior;  
    
  case 'SSGP'        
    % ============================================================
    % SSGP
    % ============================================================
    % Predictions with sparse spectral sampling approximation for GP
    % The approximation is proposed by M. Lazaro-Gredilla, J. 
    % Quinonero-Candela and A. Figueiras-Vidal in Microsoft
    % Research technical report MSR-TR-2007-152 (November 2007)
    % NOTE! This does not work at the moment.

    % First evaluate the needed covariance matrices
    % v defines that parameter is a vector
    [Phi, S] = gp_trcov(gp, x);        % n x m and nxn sparse matrices
    Sv = diag(S);
    
    m = size(Phi,2);
    
    A = eye(m,m) + Phi'*(S\Phi);
    [LA, notpositivedefinite] = chol(A,'lower');
    if notpositivedefinite
      % If not positive definite, return NaN
      g=NaN; gdata=NaN; gprior=NaN;
      return;
    end
    L = (S\Phi)/A';

    b = y'./Sv' - (y'*L)*L';
    iSPhi = S\Phi;
    
    % =================================================================
    if ~isempty(strfind(gp.infer_params,'covariance'))
      % Loop over the covariance functions
      for i=1:ncf
        i1=0;
        if ~isempty(gprior)
          i1 = length(gprior);
        end
        
        gpcf = gp.cf{i};
        
        
        % Get the gradients of the covariance matrices 
        % and gprior from gpcf_* structures
        DKff = gpcf.fh.cfg(gpcf, x);
        gprior_cf = -gpcf.fh.lpg(gpcf);

        % Evaluate the gradient with respect to lengthScale
        for i2 = 1:length(DKff)
          i1 = i1+1;
          iSDPhi = S\DKff{i2};
          
          gdata(i1) = 0.5*( sum(sum(iSDPhi.*Phi,2)) + sum(sum(iSPhi.*DKff{i2},2)) );
          gdata(i1) = gdata(i1) - 0.5*( sum(sum(L'.*(L'*DKff{i2}*Phi' + L'*Phi*DKff{i2}'),1)) );
          gdata(i1) = gdata(i1) - 0.5*(b*DKff{i2}*Phi' + b*Phi*DKff{i2}')*b';
          gprior(i1) = gprior_cf(i2);
        end

        % Set the gradients of hyperparameter
        if length(gprior_cf) > length(DKff)
          for i2=length(DKff)+1:length(gprior_cf)
            i1 = i1+1;
            gdata(i1) = 0;
            gprior(i1) = gprior_cf(i2);
          end
        end        
      end
    end
    
    % =================================================================
    % Gradient with respect to Gaussian likelihood function parameters
    if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh,'trcov')
      % Evaluate the gradient from Gaussian likelihood
      DCff = gp.lik.fh.cfg(gp.lik, x);
      gprior_lik = -gp.lik.fh.lpg(gp.lik);
      for i2 = 1:length(DCff)
        i1 = i1+1;
        gdata(i1)= -0.5*DCff{i2}.*b*b';
        gdata(i1)= gdata(i1) + 0.5*sum(1./Sv-sum(L.*L,2)).*DCff{i2};
        gprior(i1) = gprior_lik(i2);
      end
      
      % Set the gradients of hyperparameter                                
      if length(gprior_lik) > length(DCff)
        for i2=length(DCff)+1:length(gprior_lik)
          i1 = i1+1;
          gdata(i1) = 0;
          gprior(i1) = gprior_lik(i2);
        end
      end                
    end
    
    % =================================================================
    % Gradient with respect to inducing inputs
    if ~isempty(strfind(gp.infer_params, 'inducing'))
      for i=1:ncf
        i1=0;
        if ~isempty(gprior)
          i1 = length(gprior);
        end
        
        gpcf = gp.cf{i};
        
        gpcf.GPtype = gp.type;        
        [gprior_ind, DKuu, DKuf] = gpcf.fh.gind(gpcf, x, y, g_ind, gdata_ind, gprior_ind);
        
        for i2 = 1:length(DKuu)
          KfuiKuuKuu = iKuuKuf'*DKuu{i2};
          
          gdata_ind(i2) = gdata_ind(i2) - 0.5.*((2*b*DKuf{i2}'-(b*KfuiKuuKuu))*(iKuuKuf*b') + ...
                                                2.*sum(sum(L'.*(L'*DKuf{i2}'*iKuuKuf))) - sum(sum(L'.*((L'*KfuiKuuKuu)*iKuuKuf))));
          gdata_ind(i2) = gdata_ind(i2) + 0.5.*(2.*b.*sum(DKuf{i2}'.*iKuuKuf',2)'*b'- b.*sum(KfuiKuuKuu.*iKuuKuf',2)'*b');
          gdata_ind(i2) = gdata_ind(i2) + 0.5.*(2.*sum(sum(L.*L,2).*sum(DKuf{i2}'.*iKuuKuf',2)) - ...
                                                sum(sum(L.*L,2).*sum(KfuiKuuKuu.*iKuuKuf',2)));                    
        end
      end
    end
    
    g = gdata + gprior;

end
