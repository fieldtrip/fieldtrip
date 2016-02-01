function [f, energ, diagn] = scaled_hmc(f, opt, gp, x, y, z)
%SCALED_HMC  A scaled hybrid Monte Carlo sampling for latent values
%
%   Description
%    [F, ENERG, DIAG] = SCALED_HMC(F, OPT, GP, X, Y) takes the
%    current latent values F, options structure OPT, Gaussian
%    process structure GP, inputs X and outputs Y. Samples new
%    latent values and returns also energies ENERG and diagnostics
%    DIAG. The latent values are sampled from their conditional
%    posterior p(f|y,th).
%
%    The latent values are whitened with an approximate posterior
%    covariance before the sampling. This reduces the
%    autocorrelation and speeds up the mixing of the sampler. See
%    Vanhatalo and Vehtari (2007) for details on implementation.
%
%    The options structure needs to be similar to hybrid
%    Monte Carlo, HMC2. See HMC2 and HMC2_OPT for details.
%
%    OPT = SCALED_HMC() Returns default options given by HMC2_OPT and
%      repeat              - the number of HMC-chains before 
%                            returning a single sample (default 10)
%
%    OPT = SCALED_HMC(OPT) Returns default options given by
%    not yet set in OPT
%
%  See also
%    GP_MC, HMC2, HMC2_OPT
  
% Copyright (c) 2006 Aki Vehtari
% Copyright (c) 2006-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Set default options using hmc2_opt
  if nargin<=1
    if nargin==0
      f=hmc2_opt();
    else
      f=hmc2_opt(f);
    end
    if ~isfield(f,'repeat')
      f.repeat=10;
    end
    return
  end
  
% Check if an old sampler state is provided
  if isfield(opt, 'rstate')
    if ~isempty(opt.rstate)
      latent_rstate = opt.rstate;
    end
  else
    latent_rstate = sum(100*clock);
  end
  
  % Initialize variables
  [n,nin] = size(x);
  switch gp.type
    case 'FULL'
      u = [];
    case 'FIC'
      u = gp.X_u;
      Lav=[];
    case 'CS+FIC'
      u = gp.X_u;
      Labl=[];
      Lp = [];            
    case {'PIC' 'PIC_BLOCK'}
      u = gp.X_u;
      ind = gp.tr_index;
      Labl=[];
      Lp = [];
  end
  if isfield(gp.lik, 'nondiagW')
    switch gp.lik.type
      case {'LGP', 'LGPC'}
        % Do nothing
      case {'Softmax', 'Multinom'}
        [n,nout] = size(y);
      otherwise
        nout=length(gp.comp_cf);
        if isfield(gp.lik,'xtime')
          xtime=gp.lik.xtime;
          ntime = size(xtime,1);
          nl=[ntime n];
        else
          nl=repmat(n,1,length(gp.comp_cf));
        end
        nlp=length(nl); % number of latent processes
    end
  end
  n=length(y);

  J = [];
  U = [];
  iJUU = [];
  Linv=[];
  L2=[];
  H_m=[];
  test=[];
  iLaKfuic=[];
  mincut = -300;

  % Transform the latent values
  switch gp.type
    case 'FULL'
      getL(f, gp, x, y, u, z);             % Evaluate the help matrix for transformation
      w = (L2\f)';                   % Rotate f towards prior
    case 'FIC'
      getL(f, gp, x, y, u, z);          % Evaluate the help matrix for transformation
      fs = f./Lp;                    % Rotate f towards prior
      w = fs + U*((J*U'-U')*fs);     
    case {'PIC' 'PIC_BLOCK'}
      getL(f, gp, x, y, u, z);          % Evaluate the help matrix for transformation
      fs=zeros(size(f));             % Rotate f towards prior
      for i=1:length(ind)
        fs(ind{i}) = Lp{i}\f(ind{i});
      end
      w = fs + U*((J*U'-U')*fs);
    case {'CS+FIC'}
      getL(f, gp, x, y, u, z);          % Evaluate the help matrix for transformation
      fs = Lp\f;                        % Rotate f towards prior
      w = fs + U*((J*U'-U')*fs);
    otherwise 
      error('unknown type of GP\n')
  end
  
  %gradcheck(w, @f_e, @f_g, gp, x, y, u, f);
  
  % Conduct the HMC sampling for the transformed latent values
  hmc2('state',latent_rstate)
  rej = 0;
  for li=1:opt.repeat 
    [w, energ, diagn] = hmc2(@f_e, w, opt, @f_g, gp, x, y, u, z);
    w = w(end,:);
    % Find an optimal scaling during the first half of repeat
    if li<opt.repeat/2
      if diagn.rej
        opt.stepadj=max(1e-5,opt.stepadj/1.4);
      else
        opt.stepadj=min(1,opt.stepadj*1.02);
      end
    end
    rej=rej+diagn.rej/opt.repeat;
    if isfield(diagn, 'opt')
      opt=diagn.opt;
    end
  end
  w = w(end,:);
  
  % Rotate w back to the latent value space
  w=w(:);
  switch gp.type
    case 'FULL'
      f=L2*w;
    case 'FIC'
      f = Lp.*(w + U*(iJUU*w));
    case  {'PIC' 'PIC_BLOCK'}
      w2 = w + U*(iJUU*w);
      for i=1:length(ind)
        f(ind{i}) = Lp{i}*w2(ind{i});
      end
    case  {'CS+FIC'}
      w2 = w + U*(iJUU*w);
      f = Lp*w2;            
  end
  % Record the options
  opt.rstate = hmc2('state');
  diagn.opt = opt;
  diagn.rej = rej;
  diagn.lvs = opt.stepadj;

  function [g, gdata, gprior] = f_g(w, gp, x, y, u, z)
  %F_G     Evaluate gradient function for transformed GP latent values 
  %               
    
  % Force f and E to be a column vector
    w=w(:);
    
    switch gp.type
      case 'FULL'
        f = L2*w;
        f = max(f,mincut);
        gdata = - gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
        
        if ~isfield(gp,'meanf')
          b=Linv*f;
          gprior=Linv'*b;
        else
          [H_m,b_m,tmp]=mean_prep(gp,x,[]);
          M = (H_m'*b_m-f);
          b=Linv*M;
          gprior=-1*Linv'*b;
        end
        
        g = (L2'*(gdata + gprior))';
      case 'FIC'
        %        w(w<eps)=0;
        f = Lp.*(w + U*(iJUU*w));
        gdata = - gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
        f = max(f,mincut);
        gprior = f./Lav - iLaKfuic*(iLaKfuic'*f);
        g = gdata +gprior;
        g = Lp.*g;
        g = g + U*(iJUU*g);
        g = g';
      case {'PIC' 'PIC_BLOCK'}
        w2= w + U*(iJUU*w);
        for i=1:length(ind)
          f(ind{i}) = Lp{i}*w2(ind{i});
        end
        f = max(f,mincut);
        gdata = - gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
        gprior = zeros(size(gdata));
        for i=1:length(ind)
          gprior(ind{i}) = Labl{i}\f(ind{i});
        end
        gprior = gprior - iLaKfuic*(iLaKfuic'*f);
        g = gdata' + gprior';
        for i=1:length(ind)
          g(ind{i}) = g(ind{i})*Lp{i};
        end
        g = g + g*U*(iJUU);
        %g = g';
      case {'CS+FIC'}
        w2= w + U*(iJUU*w);
        f = Lp*w2;
        f = max(f,mincut);
        gdata = - gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
        gprior = zeros(size(gdata));
        gprior = ldlsolve(Labl,f);
        gprior = gprior - iLaKfuic*(iLaKfuic'*f);
        g = gdata' + gprior';
        g = g*Lp;
        g = g + g*U*(iJUU);
    end
  end

  function [e, edata, eprior] = f_e(w, gp, x, y, u, z)
  % F_E     Evaluate energy function for transformed GP latent values 
    
  % force f and E to be a column vector
    w=w(:);
    switch gp.type
      case 'FULL'
        f = L2*w;
        f = max(f,mincut);
        
        if ~isfield(gp,'meanf')
          B=Linv*f;
        else
          [H_m,b_m,B_m]=mean_prep(gp,x,[]);
          M=H_m'*b_m-f;
          
          B=Linv*M;
        end
        
        eprior=.5*sum(B.^2);
      case 'FIC' 
        f = Lp.*(w + U*(iJUU*w));
        f = max(f,mincut);                
        B = f'*iLaKfuic;  % 1 x u
        eprior = 0.5*sum(f.^2./Lav)-0.5*sum(B.^2);
      case {'PIC' 'PIC_BLOCK'}
        w2= w + U*(iJUU*w);
        for i=1:length(ind)
          f(ind{i}) = Lp{i}*w2(ind{i});
        end
        f = max(f,mincut);
        B = f'*iLaKfuic;  % 1 x u
        eprior = - 0.5*sum(B.^2);
        for i=1:length(ind)
          eprior = eprior + 0.5*f(ind{i})'/Labl{i}*f(ind{i});
        end
      case {'CS+FIC'}
        w2= w + U*(iJUU*w);
        f = Lp*w2;
        f = max(f,mincut);
        B = f'*iLaKfuic;  % 1 x u
        eprior = - 0.5*sum(B.^2);
        eprior = eprior + 0.5*f'*ldlsolve(Labl,f);
    end
    edata =  - gp.lik.fh.ll(gp.lik, y, f, z);
    e=edata + eprior;
  end

  function getL(w, gp, x, y, u, z)
  % GETL        Evaluate the transformation matrix (or matrices)
    
    if ~isfield(gp.lik, 'nondiagW')
      % Likelihoods with diagonal Hessian
      
      % Evaluate the Lambda (La) for specific model
      E = -gp.lik.fh.llg2(gp.lik, y, zeros(size(y)), 'latent', z);
      switch gp.type
        case 'FULL'
          C=gp_trcov(gp, x);
          % Evaluate a approximation for posterior variance
          % Take advantage of the matrix inversion lemma
          %        L=chol(inv(inv(C) + diag(E)))';
          
          if isfield(gp,'meanf')
            [H_m,b_m,B_m]=mean_prep(gp,x,[]);
            N = C + H_m'*B_m*H_m;
            Linv = inv(chol(N,'lower'));
          else
            Linv = inv(chol(C,'lower'));
          end
          
          L2 = C/chol(diag(1./E) + C);
          L2 = chol(C - L2*L2')';
        case 'FIC'
          [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % f x 1  vector
          K_fu = gp_cov(gp, x, u);           % f x u
          K_uu = gp_trcov(gp, u);            % u x u, noiseles covariance K_uu
          K_uu = (K_uu+K_uu')/2;             % ensure the symmetry of K_uu
          Luu = chol(K_uu)';
          % Q_ff = K_fu*inv(K_uu)*K_fu'
          % Here we need only the diag(Q_ff), which is evaluated below
          b=Luu\(K_fu');       % u x f
          Qv_ff=sum(b.^2)';
          Lav = Cv_ff-Qv_ff;   % f x 1, Vector of diagonal elements
          % Lets scale Lav to ones(f,1) so that Qff+La -> sqrt(La)*Qff*sqrt(La)+I
          % and form iLaKfu
          iLaKfu = zeros(size(K_fu));  % f x u,
          for i=1:n
            iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
          end
          c = K_uu+K_fu'*iLaKfu;
          c = (c+c')./2;         % ensure symmetry
          c = chol(c)';   % u x u,
          ic = inv(c);
          iLaKfuic = iLaKfu*ic';
          Lp = sqrt(1./(E + 1./Lav));
          b=b';
          for i=1:n
            b(i,:) = iLaKfuic(i,:).*Lp(i);
          end
          [V,S2]= eig(b'*b);
          S = sqrt(S2);
          U = b*(V/S);
          U(abs(U)<eps)=0;
          %        J = diag(sqrt(diag(S2) + 0.01^2));
          J = diag(sqrt(1-diag(S2)));   % this could be done without forming the diag matrix
          % J = diag(sqrt(2/(1+diag(S))));
          iJUU = J\U'-U';
          iJUU(abs(iJUU)<eps)=0;
        case {'PIC' 'PIC_BLOCK'}
          [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % f x 1  vector
          K_fu = gp_cov(gp, x, u);           % f x u
          K_uu = gp_trcov(gp, u);            % u x u, noiseles covariance K_uu
          K_uu = (K_uu+K_uu')/2;             % ensure the symmetry of K_uu
          Luu = chol(K_uu)';
          
          % Q_ff = K_fu*inv(K_uu)*K_fu'
          % Here we need only the diag(Q_ff), which is evaluated below
          B=Luu\(K_fu');       % u x f
          iLaKfu = zeros(size(K_fu));  % f x u
          for i=1:length(ind)
            Qbl_ff = B(:,ind{i})'*B(:,ind{i});
            [Kbl_ff, Cbl_ff] = gp_trcov(gp, x(ind{i},:));
            Labl{i} = Cbl_ff - Qbl_ff;
            iLaKfu(ind{i},:) = Labl{i}\K_fu(ind{i},:);    % Check if works by changing inv(Labl{i})!!!
          end
          % Lets scale Lav to ones(f,1) so that Qff+La -> sqrt(La)*Qff*sqrt(La)+I
          % and form iLaKfu
          A = K_uu+K_fu'*iLaKfu;
          A = (A+A')./2;            % Ensure symmetry
          
          % L = iLaKfu*inv(chol(A));
          iLaKfuic = iLaKfu*inv(chol(A));
          
          for i=1:length(ind)
            Lp{i} = chol(inv(diag(E(ind{i})) + inv(Labl{i})));
          end
          b=zeros(size(B'));
          
          for i=1:length(ind)
            b(ind{i},:) = Lp{i}*iLaKfuic(ind{i},:);
          end
          
          [V,S2]= eig(b'*b);
          S = sqrt(S2);
          U = b*(V/S);
          U(abs(U)<eps)=0;
          %        J = diag(sqrt(diag(S2) + 0.01^2));
          J = diag(sqrt(1-diag(S2)));   % this could be done without forming the diag matrix
          % J = diag(sqrt(2./(1+diag(S))));
          iJUU = J\U'-U';
          iJUU(abs(iJUU)<eps)=0;
        case 'CS+FIC'
          
          % Evaluate the FIC part of the prior covariance
          cf_orig = gp.cf;
          
          cf1 = {};
          cf2 = {};
          j = 1;
          k = 1;
          for i = 1:length(gp.cf)
            if ~isfield(gp.cf{i},'cs')
              cf1{j} = gp.cf{i};
              j = j + 1;
            else
              cf2{k} = gp.cf{i};
              k = k + 1;
            end
          end
          gp.cf = cf1;
          
          [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % n x 1  vector
          K_fu = gp_cov(gp, x, u);           % n x m
          K_uu = gp_trcov(gp, u);            % m x m, noiseles covariance K_uu
          K_uu = (K_uu+K_uu')/2;             % ensure the symmetry of K_uu
          Luu = chol(K_uu)';
          B=Luu\(K_fu');                     % m x n
          
          Qv_ff=sum(B.^2)';
          Lav = Cv_ff-Qv_ff;                 % n x 1, Vector of diagonal elements
          
          % Evaluate the CS part of the prior covariance
          gp.cf = cf2;
          K_cs = gp_trcov(gp,x);
          La = sparse(1:n,1:n,Lav,n,n) + K_cs;
          
          Labl = ldlchol(La);
          
          gp.cf = cf_orig;
          iLaKfu = ldlsolve(Labl,K_fu);
          
          % scale Lav to ones(f,1) so that Qff+La -> sqrt(La)*Qff*sqrt(La)+I
          A = K_uu+K_fu'*iLaKfu;
          A = (A+A')./2;                     % Ensure symmetry
          
          % L = iLaKfu*inv(chol(A));
          iLaKfuic = iLaKfu/chol(A);
          Lp = sparse(1:n,1:n,sqrt(1./(E + 1./diag(La))), n, n);
          
          b=zeros(size(B'));
          
          b = Lp*iLaKfuic;
          
          [V,S2]= eig(b'*b);
          S = sqrt(S2);
          U = b*(V/S);
          U(abs(U)<eps)=0;
          J = diag(sqrt(1-diag(S2)));   % this could be done without forming the diag matrix
          
          iJUU = J\U'-U';
          iJUU(abs(iJUU)<eps)=0;
      end
    else
      % Likelihoods with non-diagonal Hessian
      switch gp.lik.type
        case {'LGP' 'LGPC'}
          K = gp_trcov(gp, x);
          g2 = gp.lik.fh.llg2(gp.lik, y, zeros(size(f)), 'latent', z);
          g2sq=sqrt(g2);
          ny=sum(y);
          
          if ~isfield(gp,'meanf')
            if strcmpi(gp.lik.type,'LGPC')
              n1=gp.lik.gridn(1); n2=gp.lik.gridn(2);
              ny2=sum(reshape(y,fliplr(gp.lik.gridn)));
            end
          else
            if strcmpi(gp.lik.type,'LGPC')
              n1=gp.lik.gridn(1); n2=gp.lik.gridn(2);
              ny2=sum(reshape(y,fliplr(gp.lik.gridn)));
            end
          end
          
          if strcmpi(gp.lik.type,'LGPC')
            R=zeros(n);
            RKR=K;
            for k1=1:n1
              R((1:n2)+(k1-1)*n2,(1:n2)+(k1-1)*n2)=sqrt(ny2(k1))*(diag(g2sq((1:n2)+(k1-1)*n2))-g2((1:n2)+(k1-1)*n2)*g2sq((1:n2)+(k1-1)*n2)');
              RKR(:,(1:n2)+(k1-1)*n2)=RKR(:,(1:n2)+(k1-1)*n2)*R((1:n2)+(k1-1)*n2,(1:n2)+(k1-1)*n2);
            end
            for k1=1:n1
              RKR((1:n2)+(k1-1)*n2,:)=R((1:n2)+(k1-1)*n2,(1:n2)+(k1-1)*n2)'*RKR((1:n2)+(k1-1)*n2,:);
            end
            %RKR=R'*K*R;
            RKR(1:(n+1):end)=RKR(1:(n+1):end)+1;
            L = chol(RKR,'lower');
            
            RC=R'*K;
            
            iRC=L'\(L\RC);
            a = -R*iRC;
            a(1:(n+1):n^2)=a(1:(n+1):n^2)+1;
%             a = eye(size(g2,1)) - R*iRC;
          else
            R=-g2*g2sq'; R(1:(n+1):end)=R(1:(n+1):end)+g2sq';
            KR=bsxfun(@times,K,g2sq')-(K*g2)*g2sq';
            RKR=ny*(bsxfun(@times,g2sq,KR)-g2sq*(g2'*KR));
            RKR(1:(n+1):end)=RKR(1:(n+1):end)+1;
            L = chol(RKR,'lower');
            
            RC=bsxfun(@times,K,g2sq)-g2sq*(g2'*K);
            iRC=L'\(L\RC);
            a=-ny*(bsxfun(@times,iRC,g2sq)-g2*(g2sq'*iRC));
            a(1:(n+1):n^2)=a(1:(n+1):n^2) + 1;
            % a=eye(size(g2,1))-ny*(bsxfun(@times,iRC,g2sq)-g2*(g2sq'*iRC));
          end
          Linv=inv(chol(K,'lower'));
          L2=chol(K*a, 'lower');
          %         L2=chol(inv(inv(K)+ny*(diag(g2)-g2*g2')),'lower');
          
        case {'Softmax' 'Multinom'}
          [pi2_vec, pi2_mat] = gp.lik.fh.llg2(gp.lik, y, zeros(n,nout), 'latent', z);
          pi2 = reshape(pi2_vec,size(y));
          
          if isfield(gp, 'comp_cf')  % own covariance for each ouput component
            multicf = true;
            if length(gp.comp_cf) ~= nout
              error('SCALED_HMC2: the number of component vectors in gp.comp_cf must be the same as number of outputs.')
            end
          else
            multicf = false;
          end
          
          % Evaluate the blocks of the covariance matrix
          K = zeros(n,n,nout);
          if multicf
            for i1=1:nout
              K(:,:,i1) = gp_trcov(gp, x, gp.comp_cf{i1});
            end
          else
            for i1=1:nout
              K(:,:,i1) = gp_trcov(gp, x);
            end
          end
          
          R = repmat(1./pi2_vec,1,n).*pi2_mat;
          RE = zeros(n,n*nout);
          for i1=1:nout
            Dc=sqrt( pi2(:,i1) );
            Lc=(Dc*Dc').*K(:,:,i1);
            Lc(1:n+1:end)=Lc(1:n+1:end)+1;
            Lc=chol(Lc);
            L(:,:,i1)=Lc;
            
            Ec=Lc'\diag(Dc);
            Ec=Ec'*Ec;
            E(:,:,i1)=Ec;
            RER(:,:,i1) = R((1:n)+(i1-1)*n,:)'*Ec*R((1:n)+(i1-1)*n,:);
            RE(:,(1:n)+(i1-1)*n) = R((1:n)+(i1-1)*n,:)'*E(:,:,i1);
          end
          M=chol(sum(RER,3));
          
          % from this on conduct with full matrices
          % NOTE! Speed up should be considered in the future
          C = zeros(n*nout,n*nout);
          EE = zeros(n*nout,n*nout);
          for i1 = 1:nout
            C((1:n)+(i1-1)*n,(1:n)+(i1-1)*n) = K(:,:,i1);
            EE((1:n)+(i1-1)*n,(1:n)+(i1-1)*n) = E(:,:,i1);
          end
          Cpost = C - C*(EE-RE'*(M\(M'\RE)))*C;
          
          % Evaluate a approximation for posterior variance
          % Take advantage of the matrix inversion lemma
          %        L=chol(inv(inv(C) + diag(E)))';
          
          Linv = inv(chol(C)');
          
          %L2 = C/chol(diag(1./E) + C);
          %L2 = chol(C - L2*L2')';
          L2 = chol(Cpost)';
        otherwise
          K = zeros(sum(nl));
          if isfield(gp.lik,'xtime')
            K(1:ntime,1:ntime)=gp_trcov(gp, xtime, gp.comp_cf{1});
            K((1:n)+ntime,(1:n)+ntime) = gp_trcov(gp, x, gp.comp_cf{2});
          else
            for i1=1:nlp
              K((1:n)+(i1-1)*n,(1:n)+(i1-1)*n) = gp_trcov(gp, x, gp.comp_cf{i1});
            end
          end
          % Second derivatives of log-likelihood
          if isfield(gp.lik,'xtime')
            [llg2diag, llg2mat] = gp.lik.fh.llg2(gp.lik, y, zeros(size(f)), 'latent', z);
            Wdiag=-llg2diag; Wmat=-llg2mat;
            %W = [diag(Wdiag(1:ntime)) Wmat; Wmat' diag(Wdiag(ntime+1:end))];
          else
            %W = -gp.lik.fh.llg2(gp.lik, y, zeros(size(f)), 'latent', z);
            Wvec=-gp.lik.fh.llg2(gp.lik, y, f, 'latent',z);
            % W = [diag(Wvec(1:n,1)) diag(Wvec(1:n,2)); diag(Wvec(n+1:end,1)) diag(Wvec(n+1:end,2))]
            Wdiag=[Wvec(1:nl(1),1); Wvec(nl(1)+(1:nl(2)),2)];
          end
          KW=zeros(sum(nl));
          KW(1:nl(1),1:nl(1))=bsxfun(@times, K(1:nl(1),1:nl(1)), Wdiag(1:nl(1))');
          KW(nl(1)+(1:nl(2)),nl(1)+(1:nl(2)))=bsxfun(@times, K(nl(1)+(1:nl(2)),nl(1)+(1:nl(2))), Wdiag(nl(1)+(1:nl(2)))');
          if isfield(gp.lik,'xtime')
            KW(1:nl(1),nl(1)+(1:nl(2)))=K(1:nl(1),1:nl(1))*Wmat;
            KW(nl(1)+(1:nl(2)),1:nl(1))=K(nl(1)+(1:nl(2)),nl(1)+(1:nl(2)))*Wmat';
          else
            KW(1:nl(1),nl(1)+(1:nl(2)))=bsxfun(@times, K((1:nl(1)),(1:nl(1))), Wvec((nl(1)+1):2*n,1)');
            KW(nl(1)+(1:nl(2)),1:nl(1))=bsxfun(@times, K(nl(1)+(1:nl(2)),nl(1)+(1:nl(2))), Wvec(1:n,2)');
          end
          
          Linv=inv(chol(K, 'lower'));
          
          % B = KW + I;
          B=KW;
          B(1:nl(1)+nl(2)+1:end)=B(1:nl(1)+nl(2)+1:end) + 1;
          [l,u]=lu(B);
          L2=chol(u\(l\K),'lower');
          % L2tmp=chol(inv(inv(K)+W),'lower');
          
      end
    end      
  end
end
