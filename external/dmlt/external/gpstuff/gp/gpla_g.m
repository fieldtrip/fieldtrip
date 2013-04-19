function [g, gdata, gprior] = gpla_g(w, gp, x, y, varargin)
%GPLA_G   Evaluate gradient of Laplace approximation's marginal 
%         log posterior estimate (GPLA_E)
%
%  Description
%    G = GPLA_G(W, GP, X, Y, OPTIONS) takes a full GP parameter
%    vector W, structure GP a matrix X of input vectors and a
%    matrix Y of target vectors, and evaluates the gradient G of
%    EP's marginal log posterior estimate. Each row of X
%    corresponds to one input vector and each row of Y corresponds
%    to one target vector.
%
%    [G, GDATA, GPRIOR] = GPLA_G(W, GP, X, Y, OPTIONS) also returns
%    the data and prior contributions to the gradient.
%
%    OPTIONS is optional parameter-value pair
%      z - optional observed quantity in triplet (x_i,y_i,z_i)
%          Some likelihoods may use this. For example, in case of
%          Poisson likelihood we have z_i=E_i, that is, expected
%          value for ith case.
%  
%  See also
%    GP_SET, GP_G, GPLA_E, GPLA_PRED

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPLA_G';
  ip.addRequired('w', @(x) isvector(x) && isreal(x) && all(isfinite(x)));
  ip.addRequired('gp',@isstruct);
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.parse(w, gp, x, y, varargin{:});
  z=ip.Results.z;

  gp = gp_unpak(gp, w);       % unpak the parameters
  ncf = length(gp.cf);
  n=size(x,1);

  g = [];
  gdata = [];
  gprior = [];
  
  if isfield(gp, 'savememory') && gp.savememory
    savememory=1;
  else
    savememory=0;
  end

  % First Evaluate the data contribution to the error
  switch gp.type
    case 'FULL'
    % ============================================================
    % FULL
    % ============================================================
    
    if ~isfield(gp.lik, 'nondiagW')
      % Likelihoos with diagonal Hessian
      
      % Calculate covariance matrix and the site parameters
      K = gp_trcov(gp,x);
      if isfield(gp,'meanf')
        [H,b_m,B_m]=mean_prep(gp,x,[]);
        K=K+H'*B_m*H;
      end
      
      [e, edata, eprior, f, L, a, W, p] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
      if isnan(e)
        g=NaN; gdata=NaN; gprior=NaN;
        return;
      end
      
      if W >= 0              % This is the usual case where likelihood is log concave
        % for example, Poisson and probit
        if issparse(K)                               % use sparse matrix routines
          
          % permute
          y = y(p);
          x = x(p,:);
          K = K(p,p);
          if ~isempty(z)
            z = z(p,:);
          end
          
          sqrtW = sqrt(W);
          
          R = sqrtW*spinv(L,1)*sqrtW;
          sqrtWK = sqrtW*K;
          C = ldlsolve(L,sqrtWK);
          C2 = diag(K) - sum(sqrtWK.*C,1)';
          s2 = 0.5*C2.*gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
        else                                         % evaluate with full matrices
          sqrtW = diag(sqrt(W));
          c = L\sqrtW;
          R = c'*c;
          C2 = diag(K) - sum((c*K).^2,1)' ;
          s2 = 0.5*C2.*gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
        end
      else                         % We might end up here if the likelihood is not log-concave
        % For example Student-t likelihood.
        C = L;
        V = L*diag(W);
        R = diag(W) - V'*V;
        C2 = sum(C.^2,1)';
        s2 = 0.5*C2.*gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
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
            % If savememory option is used, just get the number of
            % hyperparameters and calculate gradients later
            np=gpcf.fh.cfg(gpcf,[],[],[],0);
          else
            DKffc = gpcf.fh.cfg(gpcf, x);
            np=length(DKffc);
          end
          gprior_cf = -gpcf.fh.lpg(gpcf);
          
          g1 = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
          for i2 = 1:np
            if savememory
              DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
            else
              DKff=DKffc{i2};
            end
            i1 = i1+1;
            if ~isfield(gp,'meanf')
              s1 = 0.5 * a'*DKff*a - 0.5*sum(sum(R.*DKff));
            else
              s1 = 0.5 * (a-K\(H'*b_m))'*DKff*(a-K\(H'*b_m)) - 0.5*sum(sum(R.*DKff));
            end
            b = DKff * g1;
            if issparse(K)
              s3 = b - K*(sqrtW*ldlsolve(L,sqrtW*b));
            else
              s3 = b - K*(R*b);
              %s3 = (1./W).*(R*b);
            end
            gdata(i1) = -(s1 + s2'*s3);
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
      if ~isempty(strfind(gp.infer_params, 'likelihood')) ...
          && ~isempty(gp.lik.fh.pak(gp.lik))
        
        gdata_lik = 0;
        lik = gp.lik;
        
        g_logPrior = -lik.fh.lpg(lik);
        if ~isempty(g_logPrior)
          
          DW_sigma = lik.fh.llg3(lik, y, f, 'latent2+param', z);
          DL_sigma = lik.fh.llg(lik, y, f, 'param', z);
          b = K * lik.fh.llg2(lik, y, f, 'latent+param', z);
          s3 = b - K*(R*b);
          nl= size(DW_sigma,2);
          
          gdata_lik = - DL_sigma - 0.5.*sum(repmat(C2,1,nl).*DW_sigma) - s2'*s3;
          
          % set the gradients into vectors that will be returned
          gdata = [gdata gdata_lik];
          gprior = [gprior g_logPrior];
          i1 = length(g_logPrior);
          i2 = length(gdata_lik);
          if i1  > i2
            gdata = [gdata zeros(1,i1-i2)];
          end
        end
      end
      
      g = gdata + gprior;
      
    else
      % Likelihoods with non-diagonal Hessian

      [n,nout]=size(y);
      if isfield(gp, 'comp_cf')  % own covariance for each ouput component
        multicf = true;
        if length(gp.comp_cf) ~= nout && nout > 1
          error('GPLA_ND_G: the number of component vectors in gp.comp_cf must be the same as number of outputs.')
        end
      else
        multicf = false;
      end
      
      % Get help parameters
      [e, edata, eprior, f, L, a, E, M] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
      if isnan(e)
        return
      end
      
      switch gp.lik.type
        
        case {'LGP', 'LGPC'}
          
          nl=n;
          nlp=length(nl); % number of latent processes
          
          if isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1
            gptmp=gp; gptmp.jitterSigma2=0;
            Ka = gp_trcov(gptmp, unique(x(:,1)));
            wtmp=gp_pak(gptmp); wtmp(1)=0; gptmp=gp_unpak(gptmp,wtmp);
            Kb = gp_trcov(gptmp, unique(x(:,2)));
            clear gptmp
            n1=size(Ka,1);
            n2=size(Kb,1);
          else
            K = gp_trcov(gp,x);
          end
          
          if isfield(gp,'meanf')
            [H,b_m,B_m]=mean_prep(gp,x,[]);
            Hb_m=H'*b_m;
            if isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1
              % only zero mean function implemented for Kronecker
              % approximation
              iKHb_m=zeros(n,1);
            else
              K=K+H'*B_m*H;
              iKHb_m=K\Hb_m;
            end
          end
          
          g2 = gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
          ny=sum(y);
          
          g3=gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
          
          if isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1
            
            [Va,Da]=eig(Ka); [Vb,Db]=eig(Kb);
            % eigenvalues of K matrix
            Dtmp=kron(diag(Da),diag(Db));
            [sDtmp,istmp]=sort(Dtmp,'descend');
            
            % Form the low-rank approximation. Exclude eigenvalues
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
            
            Lb=gp_trvar(gp,x)-sum(bsxfun(@times,Vlr.*Vlr,Dlr'),2);
            
            if isfield(gp,'meanf')
              Dt=[Dlr; diag(B_m)];
              Vt=[Vlr H'];
            else
              Dt=Dlr;
              Vt=Vlr;
            end
            
            Lbt=ny*(g2)+1./Lb;
            
            St=[diag(1./Dt)+Vt'*bsxfun(@times,1./Lb,Vt) zeros(size(Dt,1),1); ...
              zeros(1,size(Dt,1)) 1];
            Pt=[bsxfun(@times,1./Lb,Vt) sqrt(ny)*g2];
            Ptt=bsxfun(@times,1./sqrt(Lbt),Pt);
            
            StL=chol(St-Ptt'*Ptt,'lower');
            iStL=StL\(bsxfun(@times,Pt',1./Lbt'));
            
            dA=(1./Lbt+sum(iStL.*iStL)');
            iStLg3=iStL*g3;
            const1=( 0.5*ny*(sum( dA.*g3))-ny*(g3'*(g3./Lbt)+iStLg3'*iStLg3) );
            const2=(g3./Lbt)+iStL'*iStLg3;
            
            s2=const1.*g3 - 0.5*ny*dA.*g3 + ny*const2.*g3;
            
          else
            if strcmpi(gp.lik.type,'LGPC')
              n1=gp.lik.gridn(1); n2=gp.lik.gridn(2);
              ny2=sum(reshape(y,fliplr(gp.lik.gridn)));
              g2sq=sqrt(g2);
              
              R=zeros(n);
              RR=zeros(n,n2);
              for k1=1:n1
                R((1:n2)+(k1-1)*n2,(1:n2)+(k1-1)*n2)=sqrt(ny2(k1))*(diag(g2sq((1:n2)+(k1-1)*n2))-g2((1:n2)+(k1-1)*n2)*g2sq((1:n2)+(k1-1)*n2)');
                RR((1:n2)+(k1-1)*n2,:)=R((1:n2)+(k1-1)*n2,(1:n2)+(k1-1)*n2)*R((1:n2)+(k1-1)*n2,(1:n2)+(k1-1)*n2)';
              end
              KW=K;
              for k1=1:n1
                KW(:,(1:n2)+(k1-1)*n2)=KW(:,(1:n2)+(k1-1)*n2)*RR((1:n2)+(k1-1)*n2,:);
              end
              %KW=K*(R*R');
              
              KW(1:(n+1):end)=KW(1:(n+1):end)+1;
              iKW=KW\eye(n);
              A=iKW*K;
              
              s2=zeros(n,1);
              for k1=1:n1
                if ny2(k1)~=0
                  g3tmp=g3((1:n2)+(k1-1)*n2);
                  Atmp=A((1:n2)+(k1-1)*n2,(1:n2)+(k1-1)*n2);
                  for ind2=1:n2
                    g3dtmp=-g3tmp*g3tmp(ind2);
                    g3dtmp(ind2)=g3dtmp(ind2)+g3tmp(ind2);
                    s2( ind2+(k1-1)*n2 ) = -ny2(k1)*0.5*sum(diag(Atmp).*g3dtmp) ...
                      + ny2(k1)*sum(sum(Atmp.*(bsxfun(@times,g3tmp,g3dtmp'))));
                  end
                end
              end
              
            else
              KW=-(K*(sqrt(ny)*g2))*(sqrt(ny)*g2)'- bsxfun(@times, K, (-ny*g2)');
              
              KW(1:(n+1):end)=KW(1:(n+1):end)+1;
              iKW=KW\eye(n);
              A=iKW*K;
              
              const1=( 0.5*ny*(sum(A(1:(n+1):end).*g3'))-ny*sum(sum(A.*(g3*g3'))) );
              const2=sum(bsxfun(@times,A,g3));
              s2=const1.*g3 - 0.5*ny*diag(A).*g3 + ny*const2'.*g3;
            end
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
              if isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1
                
                gptmp=gp; gptmp.jitterSigma2=0;
                DKa = gpcf.fh.cfg(gptmp.cf{1}, unique(x(:,1)));
                wtmp=gp_pak(gptmp); wtmp(1)=0; gptmp=gp_unpak(gptmp,wtmp);
                DKb = gpcf.fh.cfg(gptmp.cf{1}, unique(x(:,2)));
                clear gptmp
                
                for j1=1:length(DKa)
                  [DVa{j1},DDa{j1}]=eig(DKa{j1});
                end
                for j1=1:length(DKb)
                  [DVb{j1},DDb{j1}]=eig(DKb{j1});
                end
                
                % low-rank approximation of the derivative matrix w.r.t.
                % magnitude hyperparameter, low-rank + diagonal correction
                Datmp=kron(diag(DDa{1}),diag(Db));
                [sDatmp,isatmp]=sort(Datmp,'descend');
                nlr=min([sum(sDtmp>gp.latent_opt.eig_tol) round(gp.latent_opt.eig_prct*n)]);
                Dmslr=sDatmp(1:nlr);
                Vmslr=zeros(n,nlr);
                for j1=1:nlr
                  Vmslr(:,j1)=kron(DVa{1}(:,ind(isatmp(j1),1)),Vb(:,ind(isatmp(j1),2)));
                end
                % diagonal correction
                dc=gpcf.fh.cfg(gpcf, x(1,:));
                Lms=dc{1}*ones(n,1)-sum(bsxfun(@times,Vmslr.^2,Dmslr'),2);
                
                
                % low-rank approximation of the derivative matrix w.r.t.
                % lengthscale hyperparameter, low-rank1 + lowr-rank2 + diagonal correction
                %
                % first low-rank part
                Datmp=kron(diag(DDa{2}),diag(Db));
                [sDatmp,isatmp]=sort(abs(Datmp),'descend');
                nlr=min([sum(sDtmp>gp.latent_opt.eig_tol) round(gp.latent_opt.eig_prct*n)]);
                sDatmp=Datmp(isatmp);
                Dlslr1=sDatmp(1:nlr);
                Vlslr1=zeros(n,nlr);
                for j1=1:nlr
                  Vlslr1(:,j1)=kron(DVa{2}(:,ind(isatmp(j1),1)),Vb(:,ind(isatmp(j1),2)));
                end
                % second low-rank part
                Dbtmp=kron(diag(Da),diag(DDb{2}));
                [sDbtmp,isbtmp]=sort(abs(Dbtmp),'descend');
                nlr=min([sum(sDtmp>gp.latent_opt.eig_tol) round(gp.latent_opt.eig_prct*n)]);
                sDbtmp=Dbtmp(isbtmp);
                Dlslr2=sDbtmp(1:nlr);
                Vlslr2=zeros(n,nlr);
                for j1=1:nlr
                  Vlslr2(:,j1)=kron(Va(:,ind(isbtmp(j1),1)),DVb{2}(:,ind(isbtmp(j1),2)));
                end
                % diagonal correction
                Lls=dc{2}*ones(n,1)-sum(bsxfun(@times,Vlslr1.^2,Dlslr1'),2)-sum(bsxfun(@times,Vlslr2.^2,Dlslr2'),2);
                
              else
                if savememory
                  % If savememory option is used, just get the number of
                  % hyperparameters and calculate gradients later
                  np=gpcf.fh.cfg(gpcf,[],[],[],0);
                else
                  DKffc = gpcf.fh.cfg(gpcf, x);
                  np=length(DKffc);
                end
              end
              
              gprior_cf = -gpcf.fh.lpg(gpcf);
              g1 = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
              
              if isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1
                
                for i2 = 1:length(DDa)
                  i1 = i1+1;
                  
                  if ~isfield(gp,'meanf')
                    if i2==1
                      % derivative wrt magnitude hyperparameter
                      Vmsa = Vmslr'*a;
                      s1a = 0.5*(Vmsa'*(bsxfun(@times,Dmslr,Vmsa)) + a'*(Lms.*a));
                    elseif i2==2
                      % derivative wrt lengthscale hyperparameter
                      Vls1 = Vlslr1'*a;
                      Vls2 = Vlslr2'*a;
                      s1a = 0.5*( Vls1'*bsxfun(@times,Dlslr1,Vls1) + Vls2'*bsxfun(@times,Dlslr2,Vls2) + a'*(Lls.*a));
                    end
                  else
                    if i2==1
                      % derivative wrt magnitude hyperparameter
                      Vmsa = Vmslr'*(a-iKHb_m);
                      s1a = 0.5*(Vmsa'*(bsxfun(@times,Dmslr,Vmsa)) + (a-iKHb_m)'*(Lms.*(a-iKHb_m)));
                    elseif i2==2
                      % derivative wrt lengthscale hyperparameter
                      Vls1 = Vlslr1'*(a-iKHb_m);
                      Vls2 = Vlslr2'*(a-iKHb_m);
                      s1a = 0.5*( Vls1'*bsxfun(@times,Dlslr1,Vls1) + Vls2'*bsxfun(@times,Dlslr2,Vls2) + (a-iKHb_m)'*(Lls.*(a-iKHb_m)));
                    end
                  end
                  
                  % DKg2=DKff{i2}*g2;
                  if i2==1
                    DKg2 = Vmslr*bsxfun(@times,Dmslr,Vmslr'*g2) + Lms.*g2;
                  elseif i2==2
                    DKg2 =  Vlslr1*bsxfun(@times,Dlslr1,Vlslr1'*g2) + Vlslr2*bsxfun(@times,Dlslr2,Vlslr2'*g2) + Lls.*g2;
                  end
                  
                  WDKg2=ny*(g2.*DKg2-(g2*(g2'*DKg2)));
                  s1b = -0.5*(ny)*( ( - (DKg2-((WDKg2./Lbt)+(iStL'*(iStL*WDKg2)))))'*(g2) );
                  
                  if i2==1 % magnitude hyperparameter
                    
                    % low-rank
                    WDVa=ny*( bsxfun(@times,g2,Vmslr)-g2*(g2'*Vmslr) );
                    stmp=Vmslr-(bsxfun(@times,(1./Lbt),WDVa)+(iStL'*(iStL*WDVa)));
                    s1clr = 0.5*sum( (sum(bsxfun(@times,stmp,Dmslr').*Vmslr,2))' .*(-ny*g2)' );
                    
                    % diagonal correction
                    s1cdtmp = Lms - ( ny*( (g2.*Lms)./Lbt  - (g2./Lbt).*(g2'.*Lms')' ) + ...
                      sum(iStL.* (ny*( bsxfun(@times,iStL,(g2.*Lms)') - (iStL*g2)*(g2'.*Lms') )) )' );
                    s1cd=0.5*sum( s1cdtmp' .*(-ny*g2)' );
                    s1c=s1clr+s1cd;
                    DKg = Vmslr*bsxfun(@times,Dmslr,Vmslr'*g1) + Lms.*g1;
                    
                    
                  elseif i2==2 % lengthscale hyperparameter
                    
                    % low-rank 1
                    WDVa=ny*( bsxfun(@times,g2,Vlslr1)-g2*(g2'*Vlslr1) );
                    stmp=Vlslr1-(bsxfun(@times,(1./Lbt),WDVa)+(iStL'*(iStL*WDVa)));
                    s1clr1 = 0.5*sum( (sum(bsxfun(@times,stmp,Dlslr1').*Vlslr1,2))' .*(-ny*g2)' );
                    
                    % low-rank 2
                    WDVb=ny*( bsxfun(@times,g2,Vlslr2)-g2*(g2'*Vlslr2) );
                    stmp=Vlslr2-(bsxfun(@times,(1./Lbt),WDVb)+(iStL'*(iStL*WDVb)));
                    s1clr2 = 0.5*sum( (sum(bsxfun(@times,stmp,Dlslr2').*Vlslr2,2))' .*(-ny*g2)' );
                    
                    % diagonal correction
                    s1cdtmp = Lls - ( ny*( (g2.*Lls)./Lbt  - (g2./Lbt).*(g2'.*Lls')' ) + ...
                      sum(iStL.* (ny*( bsxfun(@times,iStL,(g2.*Lls)') - (iStL*g2)*(g2'.*Lls') )) )' );
                    s1cd=0.5*sum( s1cdtmp' .*(-ny*g2)' );
                    
                    s1c=s1clr1+s1clr2+s1cd;
                    
                    DKg = Vlslr1*bsxfun(@times,Dlslr1,Vlslr1'*g1) + Vlslr2*bsxfun(@times,Dlslr2,Vlslr2'*g1) + Lls.*g1;
                    
                  end
                  
                  s1=s1a+s1b+s1c;
                  
                  %DKg=DKff{i2}*g1;
                  WDKg=ny*(g2.*DKg-(g2*(g2'*DKg)));
                  s3=DKg-((WDKg./Lbt)+(iStL'*(iStL*WDKg)));
                  
                  gdata(i1) = -(s1 + s2'*s3);
                  gprior(i1) = gprior_cf(i2);
                end
                
              else
                
                for i2 = 1:np
                  i1 = i1+1;
                  if savememory
                    DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
                  else
                    DKff=DKffc{i2};
                  end
                  if ~isfield(gp,'meanf')
                    if strcmpi(gp.lik.type,'LGPC')
                      %s1 = 0.5 * a'*DKff{i2}*a - 0.5*((-iKW*(DKff{i2}*(sqrt(ny)*g2)))'*(sqrt(ny)*g2)) + 0.5*sum(sum(iKW'.*DKff{i2}).*(-ny*g2)');
                      s1 = 0.5 * a'*DKff*a - 0.5*sum(diag( R*(R'*(iKW*DKff))));
                    else
                      s1 = 0.5 * a'*DKff*a - 0.5*((-iKW*(DKff*(sqrt(ny)*g2)))'*(sqrt(ny)*g2)) + 0.5*sum(sum(iKW'.*DKff).*(-ny*g2)');
                    end
                  else
                    if strcmpi(gp.lik.type,'LGPC')
                      s1 = 0.5 * (a-iKHb_m)'*DKff*(a-iKHb_m) - 0.5*sum(diag( R*(R'*(iKW*DKff))));
                    else
                      s1 = 0.5 * (a-iKHb_m)'*DKff*(a-iKHb_m) - 0.5*((-iKW*(DKff*(sqrt(ny)*g2)))'*(sqrt(ny)*g2)) + 0.5*sum(sum(iKW'.*DKff).*(-ny*g2)');
                    end
                  end
                  %b = DKff{i2} * g1;
                  if issparse(K)
                    s3 = b - K*(sqrtW*ldlsolve(L,sqrtW*b));
                  else
                    s3=iKW*(DKff*g1);
                  end
                  
                  gdata(i1) = -(s1 + s2'*s3);
                  gprior(i1) = gprior_cf(i2);
                end
              end
              
              if isfield(gp.latent_opt, 'kron') && gp.latent_opt.kron==1
                % Set the gradients of hyperparameter
                if length(gprior_cf) > length(DKa)
                  for i2=length(DKff)+1:length(gprior_cf)
                    i1 = i1+1;
                    gdata(i1) = 0;
                    gprior(i1) = gprior_cf(i2);
                  end
                end
              else
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
          end
          
          % =================================================================
          % Gradient with respect to likelihood function parameters
          if ~isempty(strfind(gp.infer_params, 'likelihood')) ...
              && ~isempty(gp.lik.fh.pak(gp.lik))
            
            gdata_lik = 0;
            lik = gp.lik;
            
            g_logPrior = -lik.fh.lpg(lik);
            if ~isempty(g_logPrior)
              
              DW_sigma = lik.fh.llg3(lik, y, f, 'latent2+param', z);
              DL_sigma = lik.fh.llg(lik, y, f, 'param', z);
              b = K * lik.fh.llg2(lik, y, f, 'latent+param', z);
              s3 = iKW*b;
              
              gdata_lik = - DL_sigma - 0.5.*sum(sum((A.*DW_sigma))) - s2'*s3;
              
              % set the gradients into vectors that will be returned
              gdata = [gdata gdata_lik];
              gprior = [gprior g_logPrior];
              i1 = length(g_logPrior);
              i2 = length(gdata_lik);
              if i1  > i2
                gdata = [gdata zeros(1,i1-i2)];
              end
            end
          end
          
          % g = gdata + gprior;
          
          
        case {'Softmax', 'Multinom'}
          
          K = zeros(n,n,nout);
          if multicf
            for i1=1:nout
              K(:,:,i1) = gp_trcov(gp, x, gp.comp_cf{i1});
            end
          else
            Ktmp=gp_trcov(gp, x);
            for i1=1:nout
              K(:,:,i1) = Ktmp;
            end
          end
          
          % softmax
          f2=reshape(f,n,nout);
          
          llg = gp.lik.fh.llg(gp.lik, y, f2, 'latent', z);
          [pi2_vec, pi2_mat] = gp.lik.fh.llg2(gp.lik, y, f2, 'latent', z);
          % W = -diag(pi2_vec) + pi2_mat*pi2_mat', where
          % W_ij = -d^2(log(p(y|f)))/(df_i)(df_j)
          R = repmat(1./pi2_vec,1,n).*pi2_mat;
          RE = zeros(n,n*nout);
          for i1=1:nout
            RE(:,(1:n)+(i1-1)*n) = R((1:n)+(i1-1)*n,:)'*E(:,:,i1);
          end
          
          inv_iWK=zeros(n,n,nout);
          
          % Matrices for computing the derivative of determinant term w.r.t. f
          A=zeros(nout, nout, n);
          Minv=M\(M'\eye(n));
          Minv=(Minv+Minv')./2;
          for cc1=1:nout
            EMinv=RE(:,(1:n)+(cc1-1)*n)'*Minv;
            KEMinv=K(:,:,cc1)*EMinv;
            for cc2=1:nout
              if cc2>=cc1
                if cc1==cc2
                  EMtmp = - EMinv*RE(:,(1:n)+(cc2-1)*n);
                  EMtmp = EMtmp + E(:,:,cc1);
                  inv_iWK(:,:,cc1) = EMtmp;
                  A(cc1,cc1,:) = diag(K(:,:,cc1))-sum((K(:,:,cc1)*EMtmp).*K(:,:,cc1),2);
                else
                  EMtmp = - KEMinv*RE(:,(1:n)+(cc2-1)*n);
                  A(cc1,cc2,:) = -sum(EMtmp.*K(:,:,cc2),2);
                  A(cc2,cc1,:) = -sum(EMtmp.*K(:,:,cc2),2);
                end
              end
            end
          end
          
          % Derivative of determinant term w.r.t. f
          s2=zeros(n*nout,1);
          dw_mat = gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
          for cc3=1:nout
            for ii1=1:n
              % s2(i)=-0.5*trace(inv(inv(K)+W)*dW/df_i)
              s2(ii1+(cc3-1)*n) = -0.5*trace(A(:,:,ii1)*dw_mat(:,:,cc3,ii1));
            end
          end
          
          % Loop over the covariance functions
          for i=1:ncf
            DKllg=zeros(size(a));
            EDKllg=zeros(size(a));
            DKffba=zeros(n*nout,1);
            
            % check in which components the covariance function is present
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
            
            i1=0;
            if ~isempty(gprior)
              i1 = length(gprior);
            end
            
            % Gradients from covariance functions
            gpcf = gp.cf{i};
            % DKff{j} = dK(x,x)/dtheta_j
            if savememory
              % If savememory option is used, just return number of
              % hyperparameters and calculate gradients later
              np=gpcf.fh.cfg(gpcf,[],[],[],0);
            else
              DKff = gpcf.fh.cfg(gpcf, x);
              np=length(DKff);
            end
            gprior_cf = -gpcf.fh.lpg(gpcf);
            
            for i2 = 1:np
              i1 = i1+1;
              if savememory
                DKffb=gpcf.fh.cfg(gpcf,x,[],[],i2);
              else
                DKffb=DKff{i2};
              end
              
              % Derivative of explicit terms
              trace_sum_tmp=0;
              for z1=1:nout
                if do(z1)
                  DKffba((1:n)+(z1-1)*n)=DKffb*a((1:n)+(z1-1)*n);
                  trace_sum_tmp = trace_sum_tmp + sum(sum( inv_iWK(:,:,z1) .* DKffb ));
                end
              end
              % s1=0.5*f'*inv(K)*dK/dtheta_j*inv(K)*f - 0.5*trace(inv(inv(W)+K)*dK/dtheta_j)
              s1 = 0.5 * a'*DKffba - 0.5.*trace_sum_tmp;
              
              % Derivative of f w.r.t. theta
              for z1=1:nout
                if do(z1)
                  DKllg((1:n)+(z1-1)*n)=DKffb*llg((1:n)+(z1-1)*n);
                  EDKllg((1:n)+(z1-1)*n)=E(:,:,z1)*DKllg((1:n)+(z1-1)*n);
                end
              end
              s3 = EDKllg - RE'*(M\(M'\(RE*DKllg)));
              for z1=1:nout
                s3((1:n)+(z1-1)*n)=K(:,:,z1)*s3((1:n)+(z1-1)*n);
              end
              % s3=inv(I+KW)*dK/dtheta_j*d(log(p(y|f)))/df
              s3 = DKllg - s3;
              
              gdata(i1) = -(s1 + s2'*s3);
              gprior(i1) = gprior_cf(i2);
              
            end
            
            % Set the gradients of hyper-hyperparameter
            if length(gprior_cf) > np
              for i2=np+1:length(gprior_cf)
                i1 = i1+1;
                gdata(i1) = 0;
                gprior(i1) = gprior_cf(i2);
              end
            end
          end
          
          %         % =================================================================
          %         % Gradient with respect to likelihood function parameters
          %         if ~isempty(strfind(gp.infer_params, 'likelihood')) ...
          %             && ~isempty(gp.lik.fh.pak(gp.lik))
          %
          %             gdata_likelih = 0;
          %             lik = gp.lik;
          %
          %             g_logPrior = feval(lik.fh.gprior, lik);
          %             if ~isempty(g_logPrior)
          %
          %                 DW_sigma = feval(lik.fh.llg3, lik, y, f, 'latent2+hyper', z);
          %                 DL_sigma = feval(lik.fh.llg, lik, y, f, 'hyper', z);
          %                 b = K * feval(lik.fh.llg2, lik, y, f, 'latent+hyper', z);
          %                 s3 = b - K*(R*b);
          %                 nl= size(DW_sigma,2);
          %
          %                 gdata_lik = - DL_sigma - 0.5.*sum(repmat(C2,1,nl).*DW_sigma) - s2'*s3;
          %
          %                 % set the gradients into vectors that will be returned
          %                 gdata = [gdata gdata_lik];
          %                 gprior = [gprior g_logPrior];
          %                 i1 = length(g_logPrior);
          %                 i2 = length(gdata_lik);
          %                 if i1  > i2
          %                     gdata = [gdata zeros(1,i1-i2)];
          %                 end
          %             end
          %         end
          
        otherwise
          
          if isfield(gp.lik,'xtime')
            xtime=gp.lik.xtime;
            if isfield(gp.lik, 'stratificationVariables')
              ebc_ind=gp.lik.stratificationVariables;
              ux = unique(x(:,ebc_ind), 'rows');
              gp.lik.n_u = size(ux,1);
              for i1=1:size(ux,1)
                gp.lik.stratind{i1}=(x(:,ebc_ind)==ux(i1));
              end
              [xtime1, xtime2] = meshgrid(ux, xtime);
              xtime = [xtime2(:) xtime1(:)];
              if isfield(gp.lik, 'removeStratificationVariables') && gp.lik.removeStratificationVariables
                x(:,ebc_ind)=[];
              end
            end
            ntime = size(xtime,1);
            nl=[ntime n];
            
            % Second derivatives of log-likelihood
            [llg2vec, llg2mat] = gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
            % W = [diag(Wdiag(1:ntime)) Wmat; Wmat' diag(Wdiag(ntime+1:end))]
            Wdiag=-llg2vec; Wmat=-llg2mat;
          else
            nl=repmat(n,1,length(gp.comp_cf));
            
            % Second derivatives of log-likelihood
            Wvec=-gp.lik.fh.llg2(gp.lik, y, f, 'latent',z);
            % W = [diag(Wvec(1:n,1)) diag(Wvec(1:n,2)); diag(Wvec(n+1:end,1)) diag(Wvec(n+1:end,2))]
            Wdiag=[Wvec(1:nl(1),1); Wvec(nl(1)+(1:nl(2)),2)];
          end
          nlp=length(nl); % Number of latent processes
          
          % K is block-diagonal covariance matrix where blocks correspond to
          % latent processes
          K = zeros(sum(nl));
          if isfield(gp.lik,'xtime')
            K(1:ntime,1:ntime)=gp_trcov(gp, xtime, gp.comp_cf{1});
            K((1:n)+ntime,(1:n)+ntime) = gp_trcov(gp, x, gp.comp_cf{2});
          else
            for i1=1:nlp
              K((1:n)+(i1-1)*n,(1:n)+(i1-1)*n) = gp_trcov(gp, x, gp.comp_cf{i1});
            end
          end
          
          if isfield(gp,'meanf')
            [H,b_m,B_m]=mean_prep(gp,x,[]);
            Hb_m=H'*b_m;
            K=K+H'*B_m*H;
            iKHb_m=K\Hb_m;
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
          
          % B = (I + K*W)
          B=KW; B(1:(nl(1)+nl(2)+1):end)=B(1:(nl(1)+nl(2)+1):end)+1;
          
          iB=B\eye(sum(nl));
          
          % A = inv(I+K*W)*K
          A=iB*K;
          
          s2=zeros(sum(nl),1);
          
          if isfield(gp.lik,'xtime')
            A_diag=diag(A);
            A_mat=A(1:ntime,ntime+(1:n));
            for i1=1:sum(nl)
              % Third derivatives
              [dw_diag,dw_mat]=gp.lik.fh.llg3(gp.lik, y, f, 'latent', z, i1);
              % s2(i) = -0.5*trace(inv(inv(K)+W)*dW/df_i)
              s2(i1) = 0.5*(sum(A_diag.*dw_diag)+2*sum(sum(A_mat.*dw_mat)));
            end
          else
            % Third derivatives
            dw_mat = gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
            for i1=1:n
              % s2(i) = -0.5*trace(inv(inv(K)+W)*dW/df_i)
              s2(i1) = 0.5*trace(A(i1:n:(i1+n),i1:n:(i1+n))*dw_mat(:,:,1,i1));
              s2(i1+n) = 0.5*trace(A(i1:n:(i1+n),i1:n:(i1+n))*dw_mat(:,:,2,i1));
            end
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
              
              % check in which components the covariance function is present
              do = false(nlp,1);
              for z1=1:nlp
                if any(gp.comp_cf{z1}==i)
                  do(z1) = true;
                end
              end
              
              if isfield(gp.lik,'xtime')
                if ~isempty(intersect(gp.comp_cf{1},i))
                  if savememory
                    % If savememory option is used, just get the number of
                    % hyperparametrs and calculate gradients later
                    np=gpcf.fh.cfg(gpcf,[],[],[],0);
                  else
                    DKffc = gpcf.fh.cfg(gpcf, xtime);
                    np=length(DKffc);
                  end
                else
                  if savememory
                    % If savememory option is used, just get the number of
                    % hyperparametrs and calculate gradients later
                    np=gpcf.fh.cfg(gpcf,[],[],[],0);
                  else
                    DKffc = gpcf.fh.cfg(gpcf, x);
                    np=length(DKffc);
                  end
                end
              else
                if savememory
                  % If savememory option is used, just get the number of
                  % hyperparametrs and calculate gradients later
                  np=gpcf.fh.cfg(gpcf,[],[],[],0);
                else
                  DKffc = gpcf.fh.cfg(gpcf, x);
                  np=length(DKffc);
                end
              end
              gprior_cf = -gpcf.fh.lpg(gpcf);
              g1 = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
              
              WiB11=bsxfun(@times, Wdiag(1:nl(1)),iB(1:nl(1),1:nl(1)));
              WiB12=bsxfun(@times, Wdiag(1:nl(1)),iB(1:nl(1),nl(1)+(1:nl(2))));
              WiB22=bsxfun(@times, Wdiag(nl(1)+(1:nl(2))),iB(nl(1)+(1:nl(2)),nl(1)+(1:nl(2))));
              if isfield(gp.lik,'xtime')
                WiB11=WiB11 + Wmat*iB(nl(1)+(1:nl(2)),1:nl(1));
                WiB12=WiB12 + Wmat*iB(nl(1)+(1:nl(2)),nl(1)+(1:nl(2)));
                WiB22=WiB22 + Wmat'*iB(1:nl(1),nl(1)+(1:nl(2)));
              else
                WiB11=WiB11 + bsxfun(@times,Wvec(1:n,2),iB(nl(1)+(1:nl(2)),1:nl(1)));
                WiB12=WiB12 + bsxfun(@times,Wvec(1:n,2),iB(nl(1)+(1:nl(2)),nl(1)+(1:nl(2))));
                WiB22=WiB22 + bsxfun(@times,Wvec(nl(1)+(1:nl(2)),1),iB(1:nl(1),nl(1)+(1:nl(2))));
              end
              WiB=[WiB11 WiB12; WiB12' WiB22];
              % WiB=W*inv(I+KW)
              
              for i2 = 1:np
                i1 = i1+1;
                if ~isfield(gp,'meanf')
                  dKnl = zeros(sum(nl));
                  if isfield(gp.lik,'xtime')
                    if ~isempty(intersect(gp.comp_cf{1},i)) %do(indnl)
                      if savememory
                        DKff=gpcf.fh.cfg(gpcf,xtime,[],[],i2);
                      else
                        DKff=DKffc{i2};
                      end
                      dKnl(1:ntime,1:ntime) = DKff;
                      %end
                    else
                      if savememory
                        DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
                      else
                        DKff=DKffc{i2};
                      end
                      %if do(indnl)
                      dKnl(ntime+(1:n),ntime+(1:n)) = DKff;
                      %end
                    end
                  else
                    if savememory
                      DKff=gpcf.fh.cfg(gpcf,x,[],[],i2);
                    else
                      DKff=DKffc{i2};
                    end
                    for indnl=1:nlp
                      if do(indnl)
                        dKnl((1:n)+(indnl-1)*n,(1:n)+(indnl-1)*n) = DKff;
                      end
                    end
                  end
                  % s1 = 0.5*f'*inv(K)*dK/dtheta_j*inv(K)*f -
                  % 0.5*trace(inv(inv(W)+K)*dK/dtheta_j)
                  s1 = 0.5 * a'*dKnl*a - 0.5*sum(sum((dKnl.*WiB)));
                else
                  %s1 = 0.5 * (a-K\(H'*b_m))'*DKff{i2}*(a-K\(H'*b_m)) - 0.5*sum(sum(R.*DKff{i2}));
                end
                %b = DKff{i2} * g1;
                b = dKnl*g1;
                % s3 = inv(I+KW)*dK/dtheta_j*d(log(p(y|f)))/df
                s3=iB*b;
                
                gdata(i1) = -(s1 + s2'*s3);
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
          
          % =================================================================
          % Gradient with respect to likelihood function parameters
          if ~isempty(strfind(gp.infer_params, 'likelihood')) ...
              && ~isempty(gp.lik.fh.pak(gp.lik))
            
            gdata_lik = 0;
            lik = gp.lik;
            
            g_logPrior = -lik.fh.lpg(lik);
            if ~isempty(g_logPrior)
              
              DW_sigma = lik.fh.llg3(lik, y, f, 'latent2+param', z);
              DL_sigma = lik.fh.llg(lik, y, f, 'param', z);
              b = K * lik.fh.llg2(lik, y, f, 'latent+param', z);
              s3 = iB*b;
              
              gdata_lik = - DL_sigma - 0.5.*sum(sum((A.*DW_sigma))) - s2'*s3;
              
              % set the gradients into vectors that will be returned
              gdata = [gdata gdata_lik];
              gprior = [gprior g_logPrior];
              i1 = length(g_logPrior);
              i2 = length(gdata_lik);
              if i1  > i2
                gdata = [gdata zeros(1,i1-i2)];
              end
            end
          end
          
      end
      
      
      g = gdata + gprior;
    end

    case 'FIC'
      % ============================================================
      % FIC
      % ============================================================
      g_ind = zeros(1,numel(gp.X_u));
      gdata_ind = zeros(1,numel(gp.X_u));
      gprior_ind = zeros(1,numel(gp.X_u));

      u = gp.X_u;
      m = size(u,1);

      [e, edata, eprior, f, L, a, La1] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);

      K_fu = gp_cov(gp, x, u);         % f x u
      K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
      Luu = chol(K_uu);
      B=Luu'\(K_fu');       % u x f
      iKuuKuf = Luu\B;
      
      W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
      sqrtW = sqrt(W);
      
      % Components for trace( inv(inv(W) + K) * dK) )
      Lah = 1 + sqrtW.*La1.*sqrtW;
      sWKfu = (repmat(sqrtW,1,m).*K_fu);
      B3 = repmat(Lah,1,m).\sWKfu;
      A = K_uu + sWKfu'*B3; A=(A+A')/2;
      L2 = repmat(sqrtW,1,m).*(B3/chol(A));
      iLa2W = sqrtW.*(Lah.\sqrtW);
      
      LL = sum(L2.*L2,2);
      BB = sum(B.^2)';
      
      % Evaluate s2
      C1 = L2'*B'*B;
      C2 = L2'.*repmat(La1',m,1);
      
      s2t = La1 + BB;
      s2t = s2t - (La1.*iLa2W.*La1 - sum(C2.^2)' + sum(B'.*((B*(repmat(iLa2W,1,m).*B'))*B)',2)...
                   - sum(C1.^2)' + 2*La1.*iLa2W.*BB - 2*La1.*sum(L2.*C1',2));

      s2 = 0.5*s2t.*gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
      b3 = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
      
      
      % =================================================================
      % Gradient with respect to covariance function parameters
      if ~isempty(strfind(gp.infer_params, 'covariance'))
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
            % 0.5* a'*dK*a, where a = K\f
            KfuiKuuKuu = iKuuKuf'*DKuu;
            gdata(i1) = -0.5.*((2.*a'*DKuf'-(a'*KfuiKuuKuu))*(iKuuKuf*a) + (a'.*DKff')*a...
                               - (2.*a'.*sum(DKuf'.*iKuuKuf',2)'*a-a'.*sum(KfuiKuuKuu.*iKuuKuf',2)'*a) );
            
            % trace( inv(inv(W) + K) * dQ) )
            gdata(i1) = gdata(i1) - 0.5.*sum(sum(L2'.*(2.*L2'*DKuf'*iKuuKuf - L2'*KfuiKuuKuu*iKuuKuf)));
            gdata(i1) = gdata(i1) + 0.5.*sum(DKff.*iLa2W - LL.*DKff);
            gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf'.*iKuuKuf',2)) - sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
            
            b = (2*DKuf' - KfuiKuuKuu)*(iKuuKuf*b3) + DKff.*b3 - sum((2.*DKuf'- KfuiKuuKuu).*iKuuKuf',2).*b3;
            bb = sqrtW.*(Lah.\(sqrtW.*b)) - L2*(L2'*b);
            s3 = b - (La1.*bb + B'*(B*bb));
            gdata(i1) = gdata(i1) - s2'*s3;
            
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
      
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && ~isempty(gp.lik.fh.pak(gp.lik))
        gdata_lik = 0;
        lik = gp.lik;

        
        DW_sigma = lik.fh.llg3(lik, y, f, 'latent2+param', z);
        DL_sigma = lik.fh.llg(lik, y, f, 'param', z);
        DL_f_sigma = lik.fh.llg2(lik, y, f, 'latent+param', z);
%         b = La1.*DL_f_sigma + B'*(B*DL_f_sigma);            
%         bb = (iLa2W.*b - L2*(L2'*b));
%         s3 = b - (La1.*bb + B'*(B*bb));
        b = repmat(La1,1,size(DL_f_sigma,2)).*DL_f_sigma + B'*(B*DL_f_sigma);            
        bb = (repmat(iLa2W,1,size(DL_f_sigma,2)).*b - L2*(L2'*b));
        s3 = b - (repmat(La1,1,size(DL_f_sigma,2)).*bb + B'*(B*bb));

%         gdata_lik = - DL_sigma - 0.5.*sum(s2t.*DW_sigma) - s2'*s3;
        gdata_lik = - DL_sigma - 0.5.*sum(repmat(s2t,1,size(DL_f_sigma,2)).*DW_sigma) - s2'*s3;
        
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
                i1 = i1+1;
              
                % 0.5* a'*dK*a, where a = K\f
                KfuiKuuKuu = iKuuKuf'*DKuu{i2};
                gdata(i1) = gdata(i1) -0.5.*((2.*a'*DKuf{i2}'-(a'*KfuiKuuKuu))*(iKuuKuf*a) + ...
                                           - (2.*a'.*sum(DKuf{i2}'.*iKuuKuf',2)'*a-a'.*sum(KfuiKuuKuu.*iKuuKuf',2)'*a) );
              
                % trace( inv(inv(W) + K) * dQ) )
                gdata(i1) = gdata(i1) - 0.5.*(sum(sum(L2'.*(2.*L2'*DKuf{i2}'*iKuuKuf - L2'*KfuiKuuKuu*iKuuKuf))));
                gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf{i2}'.*iKuuKuf',2)) - sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
              
              
                % b2*dK*b3
                b = (2*DKuf{i2}'-KfuiKuuKuu)*(iKuuKuf*b3)  - sum((2.*DKuf{i2}'- KfuiKuuKuu).*iKuuKuf',2).*b3;
                bb = (iLa2W.*b - L2*(L2'*b));
                s3 = b - (La1.*bb + B'*(B*bb));
                gdata(i1) = gdata(i1) - s2'*s3;
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
      m = size(u,1);
      ind = gp.tr_index;

      [e, edata, eprior, f, L, a, La1] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);

      K_fu = gp_cov(gp, x, u);         % f x u
      K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu
      Luu = chol(K_uu);
      iKuuKuf = Luu\(Luu'\K_fu');
      B=Luu'\(K_fu');       % u x f

      W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
      sqrtW = sqrt(W);
      
      % Components for trace( inv(inv(W) + K) * dK) )
      B2 = (repmat(sqrtW,1,m).*K_fu);
      for i=1:length(ind)
        La2{i} = eye(size(La1{i})) + diag(sqrtW(ind{i}))*La1{i}*diag(sqrtW(ind{i}));
        LLa2{i} = chol(La2{i});
        B3(ind{i},:) = LLa2{i}\(LLa2{i}'\B2(ind{i},:));
      end
      A2 = K_uu + B2'*B3; A2=(A2+A2')/2;
      L2 = repmat(sqrtW,1,m).*B3/chol(A2);
      for i=1:length(ind)
        iLa2W{i} = diag(sqrtW(ind{i}))*(LLa2{i}\(LLa2{i}'\diag(sqrtW(ind{i}))));
      end
      
      LL = sum(L2.*L2,2);
      BB = sum(B.^2)';
      
      % Evaluate s2
      C1 = L2'*B'*B;
      s2t = BB;
      for i=1:length(ind)
        C2(:,ind{i}) = L2(ind{i},:)'*La1{i};
        s2t1(ind{i},:) = diag(La1{i}*iLa2W{i}*La1{i});
        s2t2(ind{i},:) = La1{i}*iLa2W{i}*B(:,ind{i})';
        s2t3(ind{i},:) = La1{i}*L2(ind{i},:);
        Bt(ind{i},:) = iLa2W{i}*B(:,ind{i})';
        s2t(ind{i}) = s2t(ind{i}) + diag(La1{i});
      end
      
      s2t = s2t - (s2t1 - sum(C2.^2)' + sum(B'.*((B*Bt)*B)',2)...
                   - sum(C1.^2)' + 2*sum(s2t2.*B',2) - 2*sum(s2t3.*C1',2));

      s2 = 0.5*s2t.*gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
      b3 = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);

      % =================================================================
      % Gradient with respect to covariance function parameters
      if ~isempty(strfind(gp.infer_params, 'covariance'))
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
            gdata(i1) = -0.5.*((2.*a'*DKuf'-(a'*KfuiKuuKuu))*(iKuuKuf*a) );
            gdata(i1) = gdata(i1) - 0.5.*(sum(sum(L2'.*(2.*L2'*DKuf'*iKuuKuf - L2'*KfuiKuuKuu*iKuuKuf))));

            b = (2*DKuf'-KfuiKuuKuu)*(iKuuKuf*b3);
            for kk=1:length(ind)
              if savememory
                DKff=gpcf.fh.cfg(gpcf, x(ind{kk},:),[],[],i2);
              else
                DKff=DKffc{kk}{i2};
              end
              gdata(i1) = gdata(i1) -0.5.*(a(ind{kk})'*DKff*a(ind{kk})...
                                           - (2.*a(ind{kk})'*DKuf(:,ind{kk})'*iKuuKuf(:,ind{kk})*a(ind{kk})...
                                              -a(ind{kk})'*KfuiKuuKuu(ind{kk},:)*iKuuKuf(:,ind{kk})*a(ind{kk})) );
              
              % trace( inv(inv(W) + K) * dQ) )                        
              gdata(i1) = gdata(i1) + 0.5.*(sum(sum(iLa2W{kk}.*DKff)) - sum(sum(L2(ind{kk},:)'.*(L2(ind{kk},:)'*DKff))));
              gdata(i1) = gdata(i1) + 0.5.*(2.*sum(sum(L2(ind{kk},:)'.*(L2(ind{kk},:)'*DKuf(:,ind{kk})'*iKuuKuf(:,ind{kk})))) - ...
                                            sum(sum(L2(ind{kk},:)'.*((L2(ind{kk},:)'*KfuiKuuKuu(ind{kk},:))*iKuuKuf(:,ind{kk})))));
              
              b(ind{kk}) = b(ind{kk}) + DKff*b3(ind{kk})...
                  - (2.*DKuf(:,ind{kk})'- KfuiKuuKuu(ind{kk},:))*iKuuKuf(:,ind{kk})*b3(ind{kk});
              bbt(ind{kk},:) = iLa2W{kk}*b(ind{kk});
            end
            
            % b2*dK*b3
            bb = (bbt - L2*(L2'*b));
            for kk=1:length(ind)
              s3t(ind{kk},:) = La1{kk}*bb(ind{kk});
            end
            s3 = b - (s3t + B'*(B*bb));
            gdata(i1) = gdata(i1) - s2'*s3;
            
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
      
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && ~isempty(gp.lik.fh.pak(gp.lik))
        gdata_lik = 0;
        lik = gp.lik;
        
        DW_sigma = lik.fh.llg3(lik, y, f, 'latent2+param', z);
        DL_sigma = lik.fh.llg(lik, y, f, 'param', z);
        DL_f_sigma = lik.fh.llg2(lik, y, f, 'latent+param', z);
        b = B'*(B*DL_f_sigma);
        for kk=1:length(ind)
          b(ind{kk}) = b(ind{kk}) + La1{kk}*DL_f_sigma(ind{kk});
          bbt(ind{kk},:) = iLa2W{kk}*b(ind{kk});
        end
        bb = (bbt - L2*(L2'*b));
        for kk=1:length(ind)
          s3t(ind{kk},:) = La1{kk}*bb(ind{kk});
        end
        s3 = b - (s3t + B'*(B*bb));

        gdata_lik = - DL_sigma - 0.5.*sum(s2t.*DW_sigma) - s2'*s3;

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
              
              
                KfuiKuuKuu = iKuuKuf'*DKuu{i2};
                gdata(i1) = -0.5.*((2.*a'*DKuf{i2}'-(a'*KfuiKuuKuu))*(iKuuKuf*a) );
                gdata(i1) = gdata(i1) - 0.5.*(sum(sum(L2'.*(2.*L2'*DKuf{i2}'*iKuuKuf - L2'*KfuiKuuKuu*iKuuKuf))));
              
                b = (2*DKuf{i2}'-KfuiKuuKuu)*(iKuuKuf*b3);
                for kk=1:length(ind)
                  gdata(i1) = gdata(i1) -0.5.*(- (2.*a(ind{kk})'*DKuf{i2}(:,ind{kk})'*iKuuKuf(:,ind{kk})*a(ind{kk})...
                                                -a(ind{kk})'*KfuiKuuKuu(ind{kk},:)*iKuuKuf(:,ind{kk})*a(ind{kk})) );
                
                  % trace( inv(inv(W) + K) * dQ) )                         
                  gdata(i1) = gdata(i1) + 0.5.*(2.*sum(sum(L2(ind{kk},:)'.*(L2(ind{kk},:)'*DKuf{i2}(:,ind{kk})'*iKuuKuf(:,ind{kk})))) - ...
                                              sum(sum(L2(ind{kk},:)'.*((L2(ind{kk},:)'*KfuiKuuKuu(ind{kk},:))*iKuuKuf(:,ind{kk})))));
                
                  b(ind{kk}) = b(ind{kk}) + ...
                    - (2.*DKuf{i2}(:,ind{kk})'- KfuiKuuKuu(ind{kk},:))*iKuuKuf(:,ind{kk})*b3(ind{kk});
                  bbt(ind{kk},:) = iLa2W{kk}*b(ind{kk});
                end
              
                % b2*dK*b3
                bb = (bbt - L2*(L2'*b));
                for kk=1:length(ind)
                  s3t(ind{kk},:) = La1{kk}*bb(ind{kk});
                end
                s3 = b - (s3t + B'*(B*bb));
                gdata(i1) = gdata(i1) - s2'*s3;
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
      u = gp.X_u;
      m = size(u,1);

      [e, edata, eprior, f, L, a, La1] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);

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
      gp.cf = cf_orig;
      
      W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
      
      % Find fill reducing permutation and permute all the
      % matrices
      p = analyze(La1);
      if ~isempty(z)
        z = z(p,:);
      end
      f = f(p);
      y = y(p);
      La1 = La1(p,p);
      K_fu = K_fu(p,:);
      L = L(p,:);
      x = x(p,:);
      W = W(p);
      a = a(p);
      
      Luu = chol(K_uu)';
      B=Luu\(K_fu');       % u x f
      iKuuKuf = Luu'\B;
      sW = sqrt(W);
      sqrtW = sparse(1:n,1:n,sW,n,n);
      Inn = sparse(1:n,1:n,1,n,n);
      
      % Components for trace( inv(inv(W) + K) * dK) )
      Lah = Inn + sqrtW*La1*sqrtW;
      LD2 = ldlchol(Lah);
      B2 = (repmat(sW,1,m).*K_fu);
      %B3 = Lah\B2;
      B3 = ldlsolve(LD2,B2);
      A2 = K_uu + B2'*B3; A2=(A2+A2')/2;
      L2 = repmat(sW,1,m).*B3/chol(A2);

      siLa2 = spinv(LD2,1);
      dsiLa2 = diag(siLa2);
      
      LL = sum(L2.*L2,2);
      BB = sum(B.^2)';
      
      % Evaluate s2
      C1 = L2'*B'*B;
      C2 = L2'*La1;
      C3 = repmat(sW,1,m).*ldlsolve(LD2,repmat(sW,1,m).*B');
      
      s2t = diag(La1) + BB;        
      %diag(La1*sqrtW*ldlsolve(LD2,sqrtW*La1))
      s2t = s2t - (diag(La1) - sum(La1*sqrtW.*siLa2',2)./sW - sum(C2.^2)' + sum(B'.*(B*C3*B)',2)...
                   - sum(C1.^2)' + 2*sum((La1*C3).*B',2) - 2*sum(C2'.*C1',2));
      
      s2 = 0.5*s2t.*gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
      
      b3 = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
      
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
              % 0.5* a'*dK*a, where a = K\f
              KfuiKuuKuu = iKuuKuf'*DKuu;
              gdata(i1) = -0.5.*((2.*a'*DKuf'-(a'*KfuiKuuKuu))*(iKuuKuf*a) + (a'.*DKff')*a...
                                 - (2.*a'.*sum(DKuf'.*iKuuKuf',2)'*a-a'.*sum(KfuiKuuKuu.*iKuuKuf',2)'*a) );
              
              % trace( inv(inv(W) + K) * dQ) )
              gdata(i1) = gdata(i1) - 0.5.*(sum(sum(L2'.*(2.*L2'*DKuf'*iKuuKuf - L2'*KfuiKuuKuu*iKuuKuf))));
              gdata(i1) = gdata(i1) + 0.5.*(sum(DKff.*dsiLa2.*W - LL.*DKff));
              gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf'.*iKuuKuf',2)) - sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
              
              gdata(i1) = gdata(i1) + 0.5.*sum(sum(sqrtW*ldlsolve(LD2,repmat(sW,1,m).*(2.*DKuf' - KfuiKuuKuu)).*iKuuKuf',2));
              gdata(i1) = gdata(i1) - 0.5.*sum(sW.*dsiLa2.*sW.*sum((2.*DKuf' - KfuiKuuKuu).*iKuuKuf',2) ); 
              
              % b2*dK*b3
              b = (2*DKuf'-KfuiKuuKuu)*(iKuuKuf*b3) + DKff.*b3 - sum((2.*DKuf'- KfuiKuuKuu).*iKuuKuf',2).*b3;
              bb = (sW.*ldlsolve(LD2,sW.*b) - L2*(L2'*b));
              s3 = b - (La1*bb + B'*(B*bb));
              gdata(i1) = gdata(i1) - s2'*s3;
              
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
              
              % Evaluate the gradient with respect to magnSigma
              gdata(i1) = 0.5*(sum(sW.*sum(siLa2.*(sqrtW*DKff)',2)) - sum(sum(L2.*(L2'*DKff')')) - a'*DKff*a);
              b = DKff*b3;
              bb = (sW.*ldlsolve(LD2,sW.*b) - L2*(L2'*b));
              s3 = b - (La1*bb + B'*(B*bb));
              gdata(i1) = gdata(i1) - s2'*s3;
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
      
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && ~isempty(gp.lik.fh.pak(gp.lik))
        gdata_lik = 0;
        lik = gp.lik;
        
        DW_sigma = lik.fh.llg3(lik, y, f, 'latent2+param', z);
        DL_sigma = lik.fh.llg(lik, y, f, 'param', z);
        DL_f_sigma = lik.fh.llg2(lik, y, f, 'latent+param', z);
        b = La1*DL_f_sigma + B'*(B*DL_f_sigma);            
        bb = (sW.*ldlsolve(LD2,sW.*b) - L2*(L2'*b));
        s3 = b - (La1*bb + B'*(B*bb));            
        
        gdata_lik = - DL_sigma - 0.5.*sum(s2t.*DW_sigma) - s2'*s3;
        
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
                  i1=i1+1;
                
                  % 0.5* a'*dK*a, where a = K\f
                  KfuiKuuKuu = iKuuKuf'*DKuu{i2};
                  gdata(i1) = -0.5.*((2.*a'*DKuf{i2}'-(a'*KfuiKuuKuu))*(iKuuKuf*a) ...
                                   - (2.*a'.*sum(DKuf{i2}'.*iKuuKuf',2)'*a-a'.*sum(KfuiKuuKuu.*iKuuKuf',2)'*a) );
                
                  % trace( inv(inv(W) + K) * dQ) )
                  gdata(i1) = gdata(i1) - 0.5.*(sum(sum(L2'.*(2.*L2'*DKuf{i2}'*iKuuKuf - L2'*KfuiKuuKuu*iKuuKuf))));
                  gdata(i1) = gdata(i1) + 0.5.*(2.*sum(LL.*sum(DKuf{i2}'.*iKuuKuf',2)) - sum(LL.*sum(KfuiKuuKuu.*iKuuKuf',2)));
                
                  gdata(i1) = gdata(i1) + 0.5.*sum(sum(sqrtW*ldlsolve(LD2,repmat(sW,1,size(u,1)).*(2.*DKuf{i2}' - KfuiKuuKuu)).*iKuuKuf',2));
                  gdata(i1) = gdata(i1) - 0.5.*( sum(sW.*dsiLa2.*sW.*sum((2.*DKuf{i2}' - KfuiKuuKuu).*iKuuKuf',2)) );
                
                  % b2*dK*b3
                  b = (2*DKuf{i2}'-KfuiKuuKuu)*(iKuuKuf*b3) - sum((2.*DKuf{i2}'- KfuiKuuKuu).*iKuuKuf',2).*b3;
                  bb = (sW.*ldlsolve(LD2,sW.*b) - L2*(L2'*b));
                  s3 = b - (La1*bb + B'*(B*bb));
                  gdata(i1) = gdata(i1) - s2'*s3;
                end
              end
            end
          end
        end
      end

      g = gdata + gprior;
      
    case {'DTC', 'VAR', 'SOR'}
      % ============================================================
      % DTC, VAR, SOR
      % ============================================================
      g_ind = zeros(1,numel(gp.X_u));
      gdata_ind = zeros(1,numel(gp.X_u));
      gprior_ind = zeros(1,numel(gp.X_u));
      
      u = gp.X_u;
      m = size(u,1);
      
      [e, edata, eprior, f, L, a, La2] = gpla_e(gp_pak(gp), gp, x, y, 'z', z);
      
      K_fu = gp_cov(gp, x, u);         % f x u
      K_uu = gp_trcov(gp, u);          % u x u, noiseles covariance K_uu
      Luu = chol(K_uu);
      B=Luu'\(K_fu');       % u x f
      iKuuKuf = Luu\B;
      
      W = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
      sqrtW = sqrt(W);
      
      % Components for trace( inv(inv(W) + K) * dK) )
      sWKfu = bsxfun(@times, sqrtW, K_fu);
      % L = chol(K_uu + sWKfu'*sWKfu)
      L2 = repmat(sqrtW,1,m).*(sWKfu/L);
      
      LL = sum(L2.*L2,2);
      BB = sum(B.^2)';
            
      % Evaluate s2
      C1 = L2'*B'*B;
      C2 = zeros(size(L2'));
      
      s2t = BB;
      s2t = s2t - (sum(C2.^2)' + sum(B'.*((B*(repmat(W,1,m).*B'))*B)',2)...
        - sum(C1.^2)');
      
      g3 = gp.lik.fh.llg3(gp.lik, y, f, 'latent', z);
      s2 = 0.5*s2t.*g3;
      b3 = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
      
      
      % =================================================================
      % Gradient with respect to covariance function parameters
      if ~isempty(strfind(gp.infer_params, 'covariance'))
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
            %             DKffc = gpcf.fh.cfg(gpcf, x, [], 1);
            DKuuc = gpcf.fh.cfg(gpcf, u);
            DKufc = gpcf.fh.cfg(gpcf, u, x);
            if isequal(gp.type, 'VAR')
              DKffc=gpcf.fh.cfg(gpcf,x,[],1);
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
            % 0.5* a'*dK*a, where a = K\f
            KfuiKuuDKuu = iKuuKuf'*DKuu;
            gdata(i1) = -0.5.*((2.*a'*DKuf'-(a'*KfuiKuuDKuu))*(iKuuKuf*a));
            
            % trace( inv(inv(W) + K) * dQ) )
            gdata(i1) = gdata(i1) - 0.5.*sum(sum(L2'.*(2.*L2'*DKuf'*iKuuKuf - L2'*KfuiKuuDKuu*iKuuKuf)));
            gdata(i1) = gdata(i1) + 0.5.*(2.*sum(W.*sum(DKuf'.*iKuuKuf',2)) - sum(W.*sum(KfuiKuuDKuu.*iKuuKuf',2)));
            
            b = (2*DKuf' - KfuiKuuDKuu)*(iKuuKuf*b3);
            bb = sqrtW.*sqrtW.*b - L2*(L2'*b);
            s3 = b - (B'*(B*bb));
            gdata(i1) = gdata(i1) - s2'*s3;

            if isequal(gp.type, 'VAR')
              % Derivative of tr(diag(K-Q).*W) (Titsias, 2009) wrt th            
              % Here we have to also split in explicit and implicit
              % derivatives
              if savememory
                DKff=gpcf.fh.cfg(gpcf,x,[],1,i2);
              else
                DKff=DKffc{i2};
              end              
              gdata(i1) = gdata(i1) + 0.5*sum(DKff.*W) - 0.5.*(2.*sum(W.*sum(DKuf'.*iKuuKuf',2)) - sum(W.*sum(KfuiKuuDKuu.*iKuuKuf',2)));
              % La2 = diag(K - Q);
              gdata(i1) = gdata(i1) - 0.5*(La2.*g3)'*s3;
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
      
      if ~isempty(strfind(gp.infer_params, 'likelihood')) && ~isempty(gp.lik.fh.pak(gp.lik))
        gdata_lik = 0;
        lik = gp.lik;
        
        
        DW_sigma = lik.fh.llg3(lik, y, f, 'latent2+param', z);
        DL_sigma = lik.fh.llg(lik, y, f, 'param', z);
        DL_f_sigma = lik.fh.llg2(lik, y, f, 'latent+param', z);
        %         b = La1.*DL_f_sigma + B'*(B*DL_f_sigma);
        %         bb = (iLa2W.*b - L2*(L2'*b));
        %         s3 = b - (La1.*bb + B'*(B*bb));
        
        % b = Kfu*inv(Kuu)*Kfu'*dW/dth
        b = B'*(B*DL_f_sigma);
        % bb = Kfu*inv(Kuu+Kfu'*W*Kfu)*Kfu'*W
        bb = repmat(1./W,1,size(DL_f_sigma,2)).*(L2*(L2'*b));
        % s3 = df/dth
        s3 = b - bb;        
        %bb = (repmat(W,1,size(DL_f_sigma,2)).*b - L2*(L2'*b));
        %s3 = b + B'*(B*bb);
        
        %         gdata_lik = - DL_sigma - 0.5.*sum(s2t.*DW_sigma) - s2'*s3;
        gdata_lik = - DL_sigma - 0.5.*sum(repmat(s2t,1,size(DL_f_sigma,2)).*DW_sigma) - s2'*s3;
        
        if isequal(gp.type, 'VAR')
          % Derivative of the trace term tr(diag(K-Q).*W)
          % Explicit dW/dth
          % La2 = diag(K-Q)
          gdata_lik = gdata_lik - 0.5*sum(La2.*DW_sigma);
          % Implicit dW/df*df/dth
          gdata_lik = gdata_lik + 0.5.*(La2.*g3)'*s3;
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
                i1 = i1+1;
                
                % 0.5* a'*dK*a, where a = K\f
                KfuiKuuDKuu = iKuuKuf'*DKuu{i2};
                gdata(i1) = -0.5.*((2.*a'*DKuf{i2}'-(a'*KfuiKuuDKuu))*(iKuuKuf*a));
                  
                % trace( inv(inv(W) + K) * dQ) )
                gdata(i1) = gdata(i1) - 0.5.*sum(sum(L2'.*(2.*L2'*DKuf{i2}'*iKuuKuf - L2'*KfuiKuuDKuu*iKuuKuf)));
                tt=0.5.*(2.*sum(W.*sum(DKuf{i2}'.*iKuuKuf',2)) - sum(W.*sum(KfuiKuuDKuu.*iKuuKuf',2)));
                gdata(i1) = gdata(i1) + tt;
                  
                b = (2*DKuf{i2}' - KfuiKuuDKuu)*(iKuuKuf*b3);
                bb = sqrtW.*sqrtW.*b - L2*(L2'*b);
                s3 = b - (B'*(B*bb));
                gdata(i1) = gdata(i1) - s2'*s3;
                
                if isequal(gp.type, 'VAR')
                  % Derivative of 0.5.*tr(diag(K-Q).*W) (Titsias, 2009) wrt X_u
                  % Here we have to also split in explicit and implicit
                  % derivatives
                  gdata(i1) = gdata(i1) + 0 - tt;
                  % La2 = diag(K - Q);
                  gdata(i1) = gdata(i1) - 0.5*(La2.*g3)'*s3;                             
                end
              end
            end
          end
        end
      end
      
      g = gdata + gprior;
      
  end
  
  assert(isreal(gdata))
  assert(isreal(gprior))
end
