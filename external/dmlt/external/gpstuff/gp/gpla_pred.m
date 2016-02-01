function [Eft, Varft, lpyt, Eyt, Varyt] = gpla_pred(gp, x, y, varargin)
%GPLA_PRED  Predictions with Gaussian Process Laplace approximation
%
%  Description
%    [EFT, VARFT] = GPLA_PRED(GP, X, Y, XT, OPTIONS)
%    takes a GP structure together with matrix X of training
%    inputs and vector Y of training targets, and evaluates the
%    predictive distribution at test inputs XT. Returns a posterior
%    mean EFT and variance VARFT of latent variables and the
%    posterior predictive mean EYT and variance VARYT.
%
%    [EFT, VARFT, LPYT] = GPLA_PRED(GP, X, Y, XT, 'yt', YT, OPTIONS)
%    returns also logarithm of the predictive density LPYT of the
%    observations YT at test input locations XT. This can be used
%    for example in the cross-validation. Here Y has to be a vector.
% 
%    [EFT, VARFT, LPYT, EYT, VARYT] = GPLA_PRED(GP, X, Y, XT, OPTIONS)
%    returns also the posterior predictive mean EYT and variance VARYT.
%
%    [EF, VARF, LPY, EY, VARY] = GPLA_PRED(GP, X, Y, OPTIONS)
%    evaluates the predictive distribution at training inputs X
%    and logarithm of the predictive density LPYT of the training
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
% Copyright (c) 2012 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPLA_PRED';
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
      if ~isfield(gp.lik, 'nondiagW')
        % Likelihoods with diagonal Hessian
        [e, edata, eprior, f, L, a, W, p] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
        
        ntest=size(xt,1);
        % notice the order xt,x to avoid transpose later
        K_nf = gp_cov(gp,xt,x,predcf);
        if isfield(gp,'meanf')
          [H,b_m,B_m,Hs]=mean_prep(gp,x,xt);
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
            Eft=Eft + K_nf*(K\H'*b_m);
          end
        end
        
        if nargout > 1
          % Evaluate the variance
          kstarstar = gp_trvar(gp,xt,predcf);
          if isfield(gp,'meanf')
            kstarstar= kstarstar + diag(Hs'*B_m*Hs);
          end
          if W >= 0
            % This is the usual case where likelihood is log concave
            % for example, Poisson and probit
            if issparse(K_nf) && issparse(L)
              % If compact support covariance functions are used
              % the covariance matrix will be sparse
              sqrtW = sqrt(W);
              sqrtWKfn = sqrtW*K_nf(:,p)';
              V = ldlsolve(L,sqrtWKfn);
              Varft = kstarstar - sum(sqrtWKfn.*V,1)';
            else
              W = diag(W);
              V = L\(sqrt(W)*K_nf');
              Varft = kstarstar - sum(V'.*V',2);
            end
          else
            % We may end up here if the likelihood is not log concace
            % For example Student-t likelihood
            V = L*diag(W);
            R = diag(W) - V'*V;
            Varft = kstarstar - sum(K_nf.*(R*K_nf')',2);
          end
        end
      else
        % Likelihoods with non-diagonal Hessian
        
        [tn,nout]=size(y);
        [e, edata, eprior, f, L, a, E, M] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
        
        switch gp.lik.type
          
          case {'LGP', 'LGPC'}
            
            W=-gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
            
            ntest=size(xt,1);
            nl=tn;
            nlt=ntest;
            nlp=length(nl); % number of latent processes
            
            if isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1
              
              gptmp=gp; gptmp.jitterSigma2=0;
              Ka = gp_trcov(gptmp, unique(x(:,1)));
              wtmp=gp_pak(gptmp); wtmp(1)=0; gptmp=gp_unpak(gptmp,wtmp);
              Kb = gp_trcov(gptmp, unique(x(:,2)));
              clear gptmp
              n1=size(Ka,1);
              n2=size(Kb,1);
              
              [Va,Da]=eig(Ka); [Vb,Db]=eig(Kb);
              
              % eigenvalues of K matrix
              Dtmp=kron(diag(Da),diag(Db));
              [sDtmp,istmp]=sort(Dtmp,'descend');
              
              n = size(y,1);
              % Form the low-rank approximation.  Exclude eigenvalues
              % smaller than gp.latent_opt.eig_tol or take
              % gp.latent_opt.eig_prct*n eigenvalues at most.
              nlr=min([sum(sDtmp>gp.latent_opt.eig_tol) round(gp.latent_opt.eig_prct*n)]);
              sDtmp=sDtmp+gp.jitterSigma2;
              
              itmp1=meshgrid(1:n1,1:n2);
              itmp2=meshgrid(1:n2,1:n1)';
              ind=[itmp1(:) itmp2(:)];
              
              % included eigenvalues
              Dlr=sDtmp(1:nlr);
              % included eigenvectors
              Vlr=zeros(n,nlr);
              for i1=1:nlr
                Vlr(:,i1)=kron(Va(:,ind(istmp(i1),1)),Vb(:,ind(istmp(i1),2)));
              end
            else
              K_nf = gp_cov(gp,xt,x,predcf);
              K = gp_trcov(gp, x);
            end
            
            if isfield(gp,'meanf')
              [H,b_m,B_m Hs]=mean_prep(gp,x,xt);
              if ~(isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1)
                K_nf=K_nf + Hs'*B_m*H;
                %K = gp_trcov(gp, x);
                K = K+H'*B_m*H;
              end
            end
            
            % Evaluate the mean
            if isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1
              Eft=f;
            else
              deriv = feval(gp.lik.fh.llg, gp.lik, y, f, 'latent', z);
              Eft = K_nf*deriv;
              if isfield(gp,'meanf')
                Eft=Eft + K_nf*(K\H'*b_m);
                %Eft=Eft + K_nf*(K\Hs'*b_m);
              end
            end
            
            if nargout > 1
              
              if isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1
                
                Lb=gp_trvar(gp,x)-sum(bsxfun(@times,Vlr.*Vlr,Dlr'),2);
                if isfield(gp,'meanf')
                  Dt=[Dlr; diag(B_m)];
                  Vt=[Vlr H'];
                else
                  Dt=Dlr;
                  Vt=Vlr;
                end
                
                g2 = feval(gp.lik.fh.llg2, gp.lik, y, f, 'latent', z);
                Lbt=sum(y)*(g2)+1./Lb;
                
                St=[diag(1./Dt)+Vt'*bsxfun(@times,1./Lb,Vt) zeros(size(Dt,1),1); ...
                  zeros(1,size(Dt,1)) 1];
                Pt=[bsxfun(@times,1./Lb,Vt) sqrt(sum(y))*g2];
                Ptt=bsxfun(@times,1./sqrt(Lbt),Pt);
                
                StL=chol(St-Ptt'*Ptt,'lower');
                iStL=StL\(bsxfun(@times,Pt',1./Lbt'));
                
                Covfd=1./Lbt;
                Covfu=iStL;
                Covf{1}=Covfd; Covf{2}=Covfu;
              else
                
                % Evaluate the variance
                if isempty(predcf)
                  kstarstarfull = gp_trcov(gp,xt);
                else
                  kstarstarfull = gp_trcov(gp,xt,predcf);
                end
                
                if isfield(gp,'meanf')
                  kstarstarfull = kstarstarfull + Hs'*B_m*Hs;
                end
                
                if strcmpi(gp.lik.type,'LGPC')
                  g2 = feval(gp.lik.fh.llg2, gp.lik, y, f, 'latent', z);
                  g2sq=sqrt(g2);
                  n1=gp.lik.gridn(1); n2=gp.lik.gridn(2);
                  ny2=sum(reshape(y,fliplr(gp.lik.gridn)));
                  
                  R=zeros(tn);
                  for k1=1:n1
                    R((1:n2)+(k1-1)*n2,(1:n2)+(k1-1)*n2)=sqrt(ny2(k1))*(diag(g2sq((1:n2)+(k1-1)*n2))-g2((1:n2)+(k1-1)*n2)*g2sq((1:n2)+(k1-1)*n2)');
                    %RKR(:,(1:n2)+(k1-1)*n2)=RKR(:,(1:n2)+(k1-1)*n2)*R((1:n2)+(k1-1)*n2,(1:n2)+(k1-1)*n2);
                  end
                  KR=K*R;
                  RKR=R'*KR;
                  RKR(1:(size(K,1)+1):end)=RKR(1:(size(K,1)+1):end)+1;
                  [L,notpositivedefinite] = chol(RKR,'lower');
                  K_nfR=K_nf*R;
                  Ltmp=L\K_nfR';
                  Covf=kstarstarfull-(Ltmp'*Ltmp);
                else
                  g2 = feval(gp.lik.fh.llg2, gp.lik, y, f, 'latent', z);
                  g2sq = sqrt(g2);
                  ny=sum(y);
                  
                  KR=bsxfun(@times,K,g2sq')-(K*g2)*g2sq';
                  RKR=ny*(bsxfun(@times,g2sq,KR)-g2sq*(g2'*KR));
                  RKR(1:(size(K,1)+1):end)=RKR(1:(size(K,1)+1):end)+1;
                  [L,notpositivedefinite] = chol(RKR,'lower');
                  
                  K_nfR=bsxfun(@times,K_nf,g2sq')-(K_nf*g2)*g2sq';
                  Ltmp=L\K_nfR';
                  Covf=kstarstarfull-ny*(Ltmp'*Ltmp);
                end
              end
            end
            
          case {'Softmax', 'Multinom'}
            
            if isfield(gp, 'comp_cf')  % own covariance for each ouput component
              multicf = true;
              if length(gp.comp_cf) ~= nout && nout > 1
                error('GPLA_ND_E: the number of component vectors in gp.comp_cf must be the same as number of outputs.')
              end
              if ~isempty(predcf)
                if ~iscell(predcf) || length(predcf)~=nout && nout > 1
                  error(['GPLA_ND_PRED: if own covariance for each output component is used,'...
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
            
            ntest=size(xt,1);
            
            % K_nf is 3-D covariance matrix where each slice corresponds to
            % each latent process (output)
            K_nf = zeros(ntest,tn,nout);
            if multicf
              for i1=1:nout
                K_nf(:,:,i1) = gp_cov(gp,xt,x,predcf{i1});
              end
            else
              for i1=1:nout
                K_nf(:,:,i1) = gp_cov(gp,xt,x,predcf{i1});
              end
            end
            
            nout=size(y,2);
            f2=reshape(f,tn,nout);
            
            llg_vec = gp.lik.fh.llg(gp.lik, y, f2, 'latent', z);
            llg = reshape(llg_vec,size(y));
            
            %mu_star = K_nf*reshape(a,tn,nout);
            a=reshape(a,size(y));
            for i1 = 1:nout
              %   Ef(:,i1) = K_nf(:,:,i1)*llg(:,i1);
              Eft(:,i1) = K_nf(:,:,i1)*a(:,i1);
            end
            
            if nargout > 1
              [pi2_vec, pi2_mat] = gp.lik.fh.llg2(gp.lik, y, f2, 'latent', z);
              % W = -diag(pi2_vec) + pi2_mat*pi2_mat', where
              % W_ij = -d^2(log(p(y|f)))/(df_i)(df_j)
              Covf=zeros(nout, nout, ntest);
              
              R=(repmat(1./pi2_vec,1,tn).*pi2_mat);
              for i1=1:nout
                b=E(:,:,i1)*K_nf(:,:,i1)';
                c_cav = R((1:tn)+(i1-1)*tn,:)*(M\(M'\(R((1:tn)+(i1-1)*tn,:)'*b)));
                
                for j1=1:nout
                  c=E(:,:,j1)*c_cav;
                  Covf(i1,j1,:)=sum(c.*K_nf(:,:,j1)');
                end
                
                kstarstar = gp_trvar(gp,xt,predcf{i1});
                Covf(i1,i1,:) = squeeze(Covf(i1,i1,:)) + kstarstar - sum(b.*K_nf(:,:,i1)')';
              end
            end
            
          otherwise
            
            ntest=size(xt,1);
            if isfield(gp.lik,'xtime')
              xtime=gp.lik.xtime;
              if isfield(gp.lik, 'stratificationVariables')
                ebc_ind=gp.lik.stratificationVariables;
                ux = unique([x(:,ebc_ind); xt(:,ebc_ind)], 'rows');
                gp.lik.n_u = size(ux,1);
                for i1=1:size(ux,1)
                  gp.lik.stratind{i1}=(x(:,ebc_ind)==ux(i1));
                  gp.lik.stratindt{i1}=(xt(:,ebc_ind)==ux(i1));
                end
                [xtime1, xtime2] = meshgrid(ux, xtime);
                xtime = [xtime2(:) xtime1(:)];
                if isfield(gp.lik, 'removeStratificationVariables') && gp.lik.removeStratificationVariables
                  x(:,ebc_ind)=[];
                  xt(:,ebc_ind)=[];
                end
              end
              ntime = size(xtime,1);
              nl=[ntime tn];
              nlt=[ntime ntest];
            else
              nl=repmat(tn,1, length(gp.comp_cf));
              nlt=repmat(ntest, 1, length(gp.comp_cf));
            end
            nlp=length(nl); % number of latent processes
            
            if isfield(gp.lik,'xtime')
              [llg2diag, llg2mat] = gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
              % W = [diag(Wdiag(1:ntime)) Wmat; Wmat' diag(Wdiag(ntime+1:end))]
              Wdiag=-llg2diag; Wmat=-llg2mat;
            else
              Wvec=-gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
              % W = [diag(Wvec(1:n,1)) diag(Wvec(1:n,2)); diag(Wvec(n+1:end,1)) diag(Wvec(n+1:end,2))]
              Wdiag=[Wvec(1:nl(1),1); Wvec(nl(1)+(1:nl(2)),2)];
            end
            
            % K_nf is K(x,xt) covariance matrix where blocks correspond to
            % latent processes
            K_nf = zeros(sum(nlt),sum(nl));
            if isempty(predcf)
              K_nf((1:nlt(2))+nlt(1),(1:nl(2))+nl(1)) = gp_cov(gp,xt,x, gp.comp_cf{2});
              if isfield(gp.lik, 'xtime')
                K_nf(1:nlt(1),1:nl(1)) = gp_cov(gp,xtime, xtime, gp.comp_cf{1});
              else
                K_nf(1:nlt(1),1:nl(1)) = gp_cov(gp,xt,x, gp.comp_cf{1});
              end
            else
              K_nf((1:nlt(2))+nlt(1),(1:nl(2))+nl(1)) = gp_cov(gp,xt,x, intersect(gp.comp_cf{2}, predcf));
              if isfield(gp.lik, 'xtime')
                K_nf(1:nlt(1),1:nl(1)) = gp_cov(gp,xtime, xtime, intersect(gp.comp_cf{1}, predcf));
              else
                K_nf(1:nlt(1),1:nl(1)) = gp_cov(gp,xt,x, intersect(gp.comp_cf{1}, predcf));
              end
            end
            
            if isfield(gp,'meanf')
              [H,b_m,B_m Hs]=mean_prep(gp,x,xt);
              if ~(isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1)
                K_nf=K_nf + Hs'*B_m*H;
                K = gp_trcov(gp, x);
                K = K+H'*B_m*H;
              end
            end
            
            deriv = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
            Eft = K_nf*deriv;
            if isfield(gp,'meanf')
              Eft=Eft + K_nf*(K\H'*b_m);
              %Eft=Eft + K_nf*(K\Hs'*b_m);
            end
            
            if nargout > 1
              
              % Evaluate the variance
              
              % Kss is K(X*,X*) covariance matrix between test points, where
              % each block corresponds to latent processes
              Kss = zeros(sum(nlt));
              if isempty(predcf)
                Kss((1:nlt(2))+nlt(1),(1:nlt(2))+nlt(1)) = gp_trcov(gp,xt,gp.comp_cf{2});
                if isfield(gp.lik,'xtime')
                  Kss(1:nlt(1),1:nlt(1)) = gp_trcov(gp,xtime,gp.comp_cf{1});
                else
                  Kss(1:nlt(1),1:nlt(1)) = gp_trcov(gp,xt,gp.comp_cf{1});
                end
              else
                Kss((1:nlt(2))+nlt(1),(1:nlt(2))+nlt(1)) = gp_trcov(gp,xt,intersect(gp.comp_cf{2}, predcf));
                if isfield(gp.lik,'xtime')
                  Kss(1:nlt(1),1:nlt(1)) = gp_trcov(gp,xtime,intersect(gp.comp_cf{1}, predcf));
                else
                  Kss(1:nlt(1),1:nlt(1)) = gp_trcov(gp,xt,intersect(gp.comp_cf{1}, predcf));
                end
              end
              
              if isfield(gp,'meanf')
                Kss = Kss + Hs'*B_m*Hs;
              end
              
              % iB = inv(I + W*K)
              iB=L\eye(sum(nl));
              
              iBW11=bsxfun(@times, iB(1:nl(1),1:nl(1)),Wdiag(1:nl(1))');
              iBW12=bsxfun(@times, iB(1:nl(1),nl(1)+(1:nl(2))), Wdiag(nl(1)+(1:nl(2)))');
              iBW22=bsxfun(@times, iB(nl(1)+(1:nl(2)),nl(1)+(1:nl(2))),Wdiag(nl(1)+(1:nl(2)))');
              if isfield(gp.lik,'xtime')
                iBW11=iBW11 + iB(1:nl(1),nl(1)+(1:nl(2)))*Wmat';
                iBW12=iBW12 + iB(1:nl(1),1:nl(1))*Wmat;
                iBW22=iBW22 + iB(nl(1)+(1:nl(2)),1:nl(1))*Wmat;
              else
                iBW11=iBW11 + bsxfun(@times,iB(1:nl(1),nl(1)+(1:nl(2))),Wvec(nl(1)+(1:nl(2)),1)');
                iBW12=iBW12 + bsxfun(@times,iB(1:nl(1),1:nl(1)),Wvec(1:nl(1),2)');
                iBW22=iBW22 + bsxfun(@times,iB(nl(1)+(1:nl(2)),1:nl(1)),Wvec(1:nl(1),2)');
              end
              iBW=[iBW11 iBW12; iBW12' iBW22];
              
              KiBWK=K_nf*iBW*K_nf';
              
              % Covf = K(X*,X*) - K(X,X*)*inv(I+WK)*W*K(X*,X)
              Covf=Kss-KiBWK;
            end
            
        end
        if nargout > 1
          Varft=Covf;
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

        % Set params for K_nf
        BB=Luu\(B');
        BB2=Luu\(K_nu');
        Varft = kstarstar - sum(BB2'.*(BB*(repmat(Lahat,1,m).\BB')*BB2)',2)  + sum((K_nu*(K_uu\(B'*L2))).^2, 2);
        
        % if the prediction is made for training set, evaluate Lav also for prediction points
        if ~isempty(tstind)
          LavsW = Lav.*sqrt(W);
          Varft(tstind) = Varft(tstind) - (LavsW./sqrt(Lahat)).^2 + sum((repmat(LavsW,1,m).*L2).^2, 2) ...
              - 2.*sum((repmat(LavsW,1,m).*(repmat(Lahat,1,m).\B)).*(K_uu\K_nu(tstind,:)')',2)...
              + 2.*sum((repmat(LavsW,1,m).*L2).*(L2'*B*(K_uu\K_nu(tstind,:)'))' ,2);
        end
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
        B = (repmat(sqrt(W),1,m).*K_fu);
        for i=1:length(ind)
          B2(ind{i},:) = Lahat{i}\B(ind{i},:);
        end
        A2 = K_uu + B'*B2; A2=(A2+A2)/2;
        L2 = B2/chol(A2);

        iKuuB = K_uu\B';
        KnfL2 = K_nu*(iKuuB*L2);
        Varft = zeros(length(xt),1);
        for i=1:length(ind)
          v_n = gp_cov(gp, xt(tstind{i},:), x(ind{i},:),predcf).*repmat(sqrtW(ind{i},:)',length(tstind{i}),1);              % n x u
          v_bu = K_nu(tstind{i},:)*iKuuB(:,ind{i});
          KnfLa = K_nu*(iKuuB(:,ind{i})/chol(Lahat{i}));
          KnfLa(tstind{i},:) = KnfLa(tstind{i},:) - (v_bu + v_n)/chol(Lahat{i});
          Varft = Varft + sum((KnfLa).^2,2);
          KnfL2(tstind{i},:) = KnfL2(tstind{i},:) - v_bu*L2(ind{i},:) + v_n*L2(ind{i},:);
        end
        Varft = kstarstar - (Varft - sum((KnfL2).^2,2));
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
      else
        tstind = [];
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
      
      K_fu = gp_cov(gp,x,u,predcf1);         % f x u
      K_uu = gp_trcov(gp,u,predcf1);    % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu
      K_nu=gp_cov(gp,xt,u,predcf1);

      Kcs_nf = gp_cov(gp, xt, x, predcf2);

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
        
        m = amd(Lahat);
        % Calculate the predictive variance according to the type
        % covariance functions used for making the prediction
        if ptype == 1 || ptype == 3                            
          % FIC part of the covariance
          Varft = kstarstar - sum(BB2'.*(BB*(Lahat\BB')*BB2)',2) + sum((K_nu*(K_uu\(B'*L2))).^2, 2);
          % Add Lav to Kcs_nf if the prediction is made for the training set
          if  ~isempty(tstind)
            % Non-CS covariance
            if ptype == 1         
              Kcs_nf = sparse(tstind,1:n,Lav,n2,n);                    
              % Non-CS and CS covariances
            else                  
              Kcs_nf = Kcs_nf + sparse(tstind,1:n,Lav,n2,n);
            end
            KcssW = Kcs_nf*sqrtW;                    
            Varft = Varft - sum((KcssW(:,m)/chol(Lahat(m,m))).^2,2) + sum((KcssW*L2).^2, 2) ...
                    - 2.*sum((KcssW*(Lahat\B)).*(K_uu\K_nu')',2) + 2.*sum((KcssW*L2).*(L2'*B*(K_uu\K_nu'))' ,2);
            % In case of both non-CS and CS prediction covariances add 
            % only Kcs_nf if the prediction is not done for the training set 
          elseif ptype == 3
            KcssW = Kcs_nf*sqrtW;
            Varft = Varft - sum((KcssW(:,m)/chol(Lahat(m,m))).^2,2) + sum((KcssW*L2).^2, 2) ...
                    - 2.*sum((KcssW*(Lahat\B)).*(K_uu\K_nu')',2) + 2.*sum((KcssW*L2).*(L2'*B*(K_uu\K_nu'))' ,2);
          end
          % Prediction with only CS covariance
        elseif ptype == 2
          KcssW = Kcs_nf*sqrtW;
          Varft = kstarstar - sum((KcssW(:,m)/chol(Lahat(m,m))).^2,2) + sum((KcssW*L2).^2, 2);
        end
      end
      
    case {'DTC' 'VAR' 'SOR'}
      % ============================================================
      % DTC, VAR, SOR
      % ============================================================
      % Predictions with DTC,VAR,SOR sparse approximation for GP
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
      
      % Evaluate the variance
      if nargout > 1
        % re-evaluate matrices with training components
        Kfu_tr = gp_cov(gp, x, u);
        Kuu_tr = gp_trcov(gp, u);
        Kuu_tr = (Kuu_tr+Kuu_tr')./2;
        
        W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
        B = bsxfun(@times, sqrt(W), Kfu_tr);
        
        % Components for (I + W^(1/2)*(Qff + La2)*W^(1/2))^(-1) = Lahat^(-1) - L2*L2'
        % L = chol(Kuu_tr + Kfu_tr'*diag(W)*Kfu_tr)
        L2 = B/L;
        
        % Set params for K_nf
        BB=Luu\(B');
        BB2=Luu\(K_nu');
        
        switch gp.type
          case 'SOR'
            Varft = sum((K_nu/Luu).^2,2)' - sum(BB2'.*(BB*(BB')*BB2)',2)  + sum((K_nu*(K_uu\(B'*L2))).^2, 2);
          case {'VAR' 'DTC'}
            kstarstar = gp_trvar(gp,xt,predcf);
            Varft = kstarstar - sum(BB2'.*(BB*(BB')*BB2)',2)  + sum((K_nu*(K_uu\(B'*L2))).^2, 2);
        end                 
      end
      
  end
  
  % ============================================================
  % Evaluate also the predictive mean and variance of new observation(s)
  % ============================================================
  if nargout == 3
    if isempty(yt)
      lpyt=[];
    else
      lpyt = gp.lik.fh.predy(gp.lik, Eft, Varft, yt, zt);
    end
  elseif nargout > 3
    [lpyt, Eyt, Varyt] = gp.lik.fh.predy(gp.lik, Eft, Varft, yt, zt);
  end
  
end
