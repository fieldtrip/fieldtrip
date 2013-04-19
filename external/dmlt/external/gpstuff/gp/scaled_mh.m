function [f, energ, diagn] = scaled_mh(f, opt, gp, x, y, z)
%SCALED_MH  A scaled Metropolis-Hastings sampling for latent values
%
%  Description
%    [F, ENERG, DIAG] = SCALED_MH(F, OPT, GP, X, Y) takes the
%    current latent values F, options structure OPT, Gaussian
%    process structure GP, inputs X and outputs Y. Samples new
%    latent values and returns also energies ENERG and diagnostics
%    DIAG. The latent values are sampled from their conditional
%    posterior p(f|y,th).
%
%    The latent values are whitened with the prior covariance
%    before the sampling. This reduces the autocorrelation and
%    speeds up the mixing of the sampler. See (Neal, 1993) for
%    details on implementation.
%
%    The options structure should include the following fields:
%      repeat              - the number of MH-steps before 
%                            returning a single sample (default 10)
%      sample_latent_scale - the scale for the MH-step (default 0.5)
%
%    OPT = SCALED_MH() Returns default options
%
%    OPT = SCALED_MH(OPT) Returns default options for fields not
%    yet set in OPT
%
%  See also
%    GP_MC
  
% Copyright (c) 1999,2011 Aki Vehtari
% Copyright (c) 2006-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% set default options using hmc2_opt
  if nargin<=1
    if nargin==0
      f=scaled_mh_opt();
    else
      f=scaled_mh_opt(f);
    end
    return
  end
  
  [n,nout] = size(y);
  if isfield(gp.lik, 'nondiagW')
    switch gp.lik.type
      case {'LGP', 'LGPC'}
        % Do nothing
      case {'Softmax', 'Multinom'}
        % Do nothing
      otherwise
        nout=length(gp.comp_cf);        
    end
    if isfield(gp, 'comp_cf')  % own covariance for each ouput component
      multicf = true;
      if length(gp.comp_cf) ~= nout
        error('SCALED_MH: the number of component vectors in gp.comp_cf must be the same as number of outputs.')
      end
    else
      multicf = false;
    end
  end
  f = reshape(f,n,nout);

  
  maxcut = -log(eps);
  mincut = -log(1/realmin - 1);
  lvs=opt.sample_latent_scale;
  a = max(min(f, maxcut),mincut);
  
  switch gp.type
    case {'FULL'}
      
      if ~isfield(gp.lik, 'nondiagW') || ismember(gp.lik.type, {'LGP' 'LGPC'})
        [K,C]=gp_trcov(gp, x);        
        if isfield(gp,'meanf')
          [H_m,b_m,B_m]=mean_prep(gp,x,[]);
          C = C + H_m'*B_m*H_m;
        end
        L=chol(C)';
      else
        L = zeros(n,n,nout);
        if multicf
          for i1=1:nout
            [tmp, C] = gp_trcov(gp, x, gp.comp_cf{i1});
            L(:,:,i1)=chol(C, 'lower');
          end
        else
          for i1=1:nout
            [tmp, C] = gp_trcov(gp, x);
            L(:,:,i1)=chol(C, 'lower');
          end
        end
      end
      e = -gp.lik.fh.ll(gp.lik, y, f, z);
      ft = zeros(size(y));
      
      % Adaptive control algorithm to find such a value for lvs 
      % that the rejection rate of Metropolis is optimal. 
      slrej = 0;
      for li=1:100
        for i1 =1:nout
          ft(:,i1)=sqrt(1-lvs.^2).*f(:,i1)+lvs.*L(:,:,i1)*randn(n,1);
        end
        ed = -gp.lik.fh.ll(gp.lik, y, ft, z);
        a=e-ed;
        if exp(a) > rand(1)
          f=ft;
          e=ed;
          lvs=min(1,lvs*1.1);
        else
          lvs=max(1e-8,lvs/1.05);
        end
      end
      opt.sample_latent_scale=lvs;
      % Do the actual sampling 
      for li=1:(opt.repeat)
        for i1 =1:nout
          ft(:,i1)=sqrt(1-lvs.^2).*f(:,i1)+lvs.*L(:,:,i1)*randn(n,1);
        end
        ed = -gp.lik.fh.ll(gp.lik, y, ft, z);
        a=e-ed;
        if exp(a) > rand(1)
          f=ft;
          e=ed;
        else
          slrej=slrej+1;
        end
      end
      diagn.rej = slrej/opt.repeat;
      diagn.lvs = lvs;
      diagn.opt=opt;
      energ=[];
      f = f(:)';
      
    case 'FIC'
      u = gp.X_u;
      m = size(u,1);
      % Turn the inducing vector on right direction
      if size(u,2) ~= size(x,2)
        u=u';
      end
      % Calculate some help matrices
      [Kv_ff, Cv_ff] = gp_trvar(gp, x);
      K_fu = gp_cov(gp, x, u);
      K_uu = gp_trcov(gp, u);
      Luu = chol(K_uu)';

      % Evaluate the Lambda (La) 
      B=Luu\(K_fu');
      Qv_ff=sum(B.^2)';
      Lav = Cv_ff-Qv_ff;
      sLav = sqrt(Lav);
      
      n=length(y);
      e = -gp.lik.fh.ll(gp.lik, y, f, z);

      % Adaptive control algorithm to find such a value for lvs 
      % so that the rejection rate of Metropolis is optimal. 
      slrej = 0;
      for li=1:100
        ft=sqrt(1-lvs.^2).*f + lvs.*(sLav.*randn(n,1) + B'*randn(m,1));
        ed = -gp.lik.fh.ll(gp.lik, y, ft, z);
        a=e-ed;
        if exp(a) > rand(1)
          f=ft;
          e=ed;
          lvs=min(1,lvs*1.1);
        else
          lvs=max(1e-8,lvs/1.05);
        end
      end
      opt.sample_latent_scale=lvs;
      % Do the actual sampling 
      for li=1:(opt.repeat)
        ft=sqrt(1-lvs.^2).*f + lvs.*(sLav.*randn(n,1) + B'*randn(m,1));
        ed = -gp.lik.fh.ll(gp.lik, y, ft, z);
        a=e-ed;
        if exp(a) > rand(1)
          f=ft;
          e=ed;
        else
          slrej=slrej+1;
        end
      end
      diagn.rej = slrej/opt.repeat;
      diagn.lvs = lvs;
      diagn.opt=opt;
      energ=[];
      f = f';        
      
    case 'PIC'
      u = gp.X_u;
      m = size(u,1);
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
        La{i} = Cbl_ff - Qbl_ff;
        CLa{i} = chol(La{i})' ;
      end
      
      n=length(y);
      e = -gp.lik.fh.ll(gp.lik, y, f, z);

      % Adaptive control algorithm to find such a value for lvs 
      % so that the rejection rate of Metropolis is optimal. 
      slrej = 0;
      for li=1:100
        sampf = randn(size(f));
        for i=1:length(ind)
          sampf(ind{i},:) = CLa{i}*sampf(ind{i},:);
        end
        ft=sqrt(1-lvs.^2).*f + lvs.*(sampf + B'*randn(m,1));
        at = max(min(ft, maxcut),mincut);
        ed = -gp.lik.fh.ll(gp.lik, y, ft, z);
        a=e-ed;
        if exp(a) > rand(1)
          f=ft;
          e=ed;
          lvs=min(1,lvs*1.1);
        else
          lvs=max(1e-8,lvs/1.05);
        end
      end
      opt.sample_latent_scale=lvs;
      % Do the actual sampling 
      for li=1:(opt.repeat)
        sampf = randn(size(f));
        for i=1:length(ind)
          sampf(ind{i},:) = CLa{i}*sampf(ind{i},:);
        end
        ft=sqrt(1-lvs.^2).*f + lvs.*(sampf + B'*randn(m,1));
        at = max(min(ft, maxcut),mincut);
        ed = -gp.lik.fh.ll(gp.lik, y, ft, z);
        a=e-ed;
        if exp(a) > rand(1)
          f=ft;
          e=ed;
        else
          slrej=slrej+1;
        end
      end
      diagn.rej = slrej/opt.repeat;
      diagn.lvs = lvs;
      diagn.opt=opt;
      energ=[];
      f = f';        
      
    case 'CS+FIC'
      u = gp.X_u;
      cf_orig = gp.cf;
      ncf = length(gp.cf);
      n = size(x,1); m = size(u,1);

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
      Luu = chol(K_uu)';

      B=Luu\(K_fu');
      Qv_ff=sum(B.^2)';
      Lav = Cv_ff-Qv_ff;   % 1 x f, Vector of diagonal elements        
      gp.cf = cf2;
      K_cs = gp_trcov(gp,x);
      La = sparse(1:n,1:n,Lav,n,n) + K_cs;
      gp.cf = cf_orig;
      
      LD = ldlchol(La);
      sLa = chol(La)';
      
      n=length(y);
      e = -gp.lik.fh.ll(gp.lik, y, f, z);

      % Adaptive control algorithm to find such a value for lvs 
      % so that the rejection rate of Metropolis is optimal. 
      slrej = 0;
      for li=1:100
        ft=sqrt(1-lvs.^2).*f + lvs.*(sLa*randn(n,1) + B'*randn(m,1));
        ed = -gp.lik.fh.ll(gp.lik, y, ft, z);
        a=e-ed;
        if exp(a) > rand(1)
          f=ft;
          e=ed;
          lvs=min(1,lvs*1.1);
        else
          lvs=max(1e-8,lvs/1.05);
        end
      end
      opt.sample_latent_scale=lvs;
      % Do the actual sampling 
      for li=1:(opt.repeat)
        ft=sqrt(1-lvs.^2).*f + lvs.*(sLa*randn(n,1) + B'*randn(m,1));
        ed = -gp.lik.fh.ll(gp.lik, y, ft, z);
        a=e-ed;
        if exp(a) > rand(1)
          f=ft;
          e=ed;
        else
          slrej=slrej+1;
        end
      end
      diagn.rej = slrej/opt.repeat;
      diagn.lvs = lvs;
      diagn.opt=opt;
      energ=[];
      f = f';
      
  end
end

function opt = scaled_mh_opt(opt)
%SCALED_MH_OPT  Default options for scaled Metropolis-Hastings sampling
%
%  Description
%    OPT = SCALED_MH_OPT
%      return default options
%    OPT = SCALED_MH_OPT(OPT)
%      fill empty options with default values
%
%  The options and defaults are
%      repeat              - the number of MH-steps before 
%                            returning a single sample (default 10)
%      sample_latent_scale - the scale for the MH-step (default 0.5)

  if nargin < 1
    opt=[];
  end

  if ~isfield(opt,'repeat')
    opt.repeat=10;
  end
  if ~isfield(opt,'sample_latent_scale')
    opt.sample_latent_scale=0.5;
  end

end
