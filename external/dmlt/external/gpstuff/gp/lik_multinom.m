function lik = lik_multinom(varargin)
%LIK_MULTINOM    Create a multinom likelihood structure 
%
%  Description
%    LIK = LIK_MULTINOM creates multinom likelihood for multi-class
%    count data. The observed numbers in each class with C classes is given
%    as 1xC vector.
%
%    The likelihood is defined as follows:
%                              __ n                __ C             
%      p(y|f^1, ..., f^C, z) = || i=1 [ gamma(N+1) || c=1 p_i^c^(y_i^c)/gamma(y_i^c+1)]
%
%    where p_i^c = exp(f_i^c)/ (sum_c=1^C exp(f_i^c)) is the succes 
%    probability for class c, which is a function of the latent variable 
%    f_i^c for the corresponding class and N=sum(y) is the number of trials.
%
%  See also
%    GP_SET, LIK_*

% Copyright (c) 2010 Jaakko Riihim�ki, Pasi Jyl�nki
% Copyright (c) 2010 Aki Vehtari
% Copyright (c) 2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_MULTINOM';
  ip.addOptional('lik', [], @isstruct);
  ip.parse(varargin{:});
  lik=ip.Results.lik;

  if isempty(lik)
    init=true;
    lik.type = 'Multinom';
    lik.nondiagW = true;
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'Multinom')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end

  if init
    % Set the function handles to the subfunctions
    lik.fh.pak = @lik_multinom_pak;
    lik.fh.unpak = @lik_multinom_unpak;
    lik.fh.ll = @lik_multinom_ll;
    lik.fh.llg = @lik_multinom_llg;    
    lik.fh.llg2 = @lik_multinom_llg2;
    lik.fh.llg3 = @lik_multinom_llg3;
    lik.fh.predy = @lik_multinom_predy;
    lik.fh.invlink = @lik_multinom_invlink;
    lik.fh.recappend = @lik_multinom_recappend;
  end

end  

function [w,s] = lik_multinom_pak(lik)
%LIK_MULTINOM_PAK  Combine likelihood parameters into one vector.
%
%  Description 
%    W = LIK_MULTINOM_PAK(LIK) takes a likelihood structure LIK and
%    returns an empty verctor W. If Multinom likelihood had
%    parameters this would combine them into a single row vector
%    W (see e.g. lik_negbin). This is a mandatory subfunction used 
%    for example in energy and gradient computations.
%     
%
%  See also
%    LIK_MULTINOM_UNPAK, GP_PAK
  
  w = []; s = {};
end


function [lik, w] = lik_multinom_unpak(lik, w)
%LIK_MULTINOM_UNPAK  Extract likelihood parameters from the vector.
%
%  Description
%    W = LIK_MULTINOM_UNPAK(W, LIK) Doesn't do anything.
% 
%    If Multinom likelihood had parameters this would extracts them
%    parameters from the vector W to the LIK structure. This is a 
%    mandatory subfunction used for example in energy and gradient 
%    computations.
%     
%
%  See also
%    LIK_MULTINOM_PAK, GP_UNPAK

  lik=lik;
  w=w;
end


function ll = lik_multinom_ll(lik, y, f, z)
%LIK_MULTINOM_LL  Log likelihood
%
%  Description
%    LL = LIK_MULTINOM_LL(LIK, Y, F) takes a likelihood structure
%    LIK, class counts Y (NxC matrix), and latent values F (NxC
%    matrix). Returns the log likelihood, log p(y|f,z). This 
%    subfunction is needed when using Laplace approximation or 
%    MCMC for inference with non-Gaussian likelihoods. This 
%    subfunction is also used in information criteria 
%    (DIC, WAIC) computations.
%
%  See also
%    LIK_MULTINOM_LLG, LIK_MULTINOM_LLG3, LIK_MULTINOM_LLG2, GPLA_E
  
  f=reshape(f,size(y));
  expf = exp(f);
  p = expf ./ repmat(sum(expf,2),1,size(expf,2));
  N = sum(y,2);
  
  ll = sum(gammaln(N+1) - sum(gammaln(y+1),2) + sum(y.*log(p),2) );
  
end


function llg = lik_multinom_llg(lik, y, f, param, z)
%LIK_MULTINOM_LLG    Gradient of the log likelihood
%
%  Description
%    LLG = LIK_MULTINOM_LLG(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, class labels Y, and latent values F. Returns
%    the gradient of the log likelihood with respect to PARAM. At
%    the moment PARAM can be 'param' or 'latent'. This subfunction 
%    is needed when using Laplace approximation or MCMC for inference 
%    with non-Gaussian likelihoods.
%
%  See also
%    LIK_MULTINOM_LL, LIK_MULTINOM_LLG2, LIK_MULTINOM_LLG3, GPLA_E
  
  f=reshape(f,size(y));
  C = size(y,2);
  expf2 = exp(f);
  N=sum(y, 2);
  pi2 = (N*ones(1,C)).*expf2./(sum(expf2, 2)*ones(1,C));
  pi_vec=pi2(:);
  llg = y(:)-pi_vec;
  
end


function [pi_vec, pi_mat] = lik_multinom_llg2(lik, y, f, param, z)
%LIK_MULTINOM_LLG2  Second gradients of the log likelihood
%
%  Description        
%    LLG2 = LIK_MULTINOM_LLG2(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, class labels Y, and latent values F. Returns
%    the Hessian of the log likelihood with respect to PARAM. At
%    the moment PARAM can be only 'latent'. LLG2 is a vector with
%    diagonal elements of the Hessian matrix (off diagonals are
%    zero). This subfunction is needed when using Laplace 
%    approximation or EP for inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_MULTINOM_LL, LIK_MULTINOM_LLG, LIK_MULTINOM_LLG3, GPLA_E
  
% multinom:
  [n,nout]=size(y);
  N = sum(y,2)*ones(1,nout);
  f=reshape(f,n,nout);
  
  expf2 = exp(f);
  pi2 = expf2./(sum(expf2, 2)*ones(1,nout));
  pi_vec=pi2(:).*N(:);
  
  pi_mat=zeros(nout*n, n);
  for i1=1:nout
    pi_mat((1+(i1-1)*n):(nout*n+1):end)=pi2(:,i1).*sqrt(N(:,i1)); 
  end
  %     D = diag(pi_vec);
  %     llg2 = -D + pi_mat*pi_mat';
  
end    

function [dw_mat] = lik_multinom_llg3(lik, y, f, param, z)
%LIK_MULTINOM_LLG3  Third gradients of the log likelihood
%
%  Description
%    LLG3 = LIK_MULTINOM_LLG3(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, class labels Y, and latent values F and
%    returns the third gradients of the log likelihood with
%    respect to PARAM. At the moment PARAM can be only 'latent'. 
%    LLG3 is a vector with third gradients. This subfunction is 
%    needed when using Laplace approximation for inference with 
%    non-Gaussian likelihoods.
%
%  See also
%    LIK_MULTINOM_LL, LIK_MULTINOM_LLG, LIK_MULTINOM_LLG2, GPLA_E, GPLA_G
  
  [n,nout] = size(y);
  f2 = reshape(f,n,nout);
  
  N=sum(y, 2);
  expf2 = exp(f2);
  pi2 = expf2./(sum(expf2, 2)*ones(1,nout));
  pi_vec=pi2(:);
  
  dw_mat=zeros(nout,nout,nout,n);
  
  for cc3=1:nout
    for ii1=1:n
      
      pic=pi_vec(ii1:n:(nout*n));
      for cc1=1:nout
        for cc2=1:nout
          
          % multinom third derivatives
          cc_sum_tmp=0;
          if cc1==cc2 && cc1==cc3 && cc2==cc3
            cc_sum_tmp=cc_sum_tmp+pic(cc1);
          end
          if cc1==cc2
            cc_sum_tmp=cc_sum_tmp-pic(cc1)*pic(cc3);
          end
          if cc2==cc3
            cc_sum_tmp=cc_sum_tmp-pic(cc1)*pic(cc2);
          end
          if cc1==cc3
            cc_sum_tmp=cc_sum_tmp-pic(cc1)*pic(cc2);
          end
          cc_sum_tmp=cc_sum_tmp+2*pic(cc1)*pic(cc2)*pic(cc3);
          
          dw_mat(cc1,cc2,cc3,ii1)=cc_sum_tmp.*N(ii1);
        end
      end
    end
  end
  
  
end

function [lpy, Ey, Vary] = lik_multinom_predy(lik, Ef, Varf, yt, zt)
%LIK_MULTINOM_PREDY  Returns the predictive mean, variance and density of y
%
%  Description
%    LPY = LIK_MULTINOM_PREDY(LIK, EF, VARF YT)
%    Returns logarithm of the predictive density PY of YT, that is
%        p(yt | y) = \int p(yt | f) p(f|y) df.
%    This requires also the incedence counts YT. This subfunction 
%    is needed when computing posterior predictive distributions for 
%    future observations.
%
%    [LPY, EY, VARY] = LIK_MULTINOM_PREDY(LIK, EF, VARF, YT) takes a
%    likelihood structure LIK, posterior mean EF and posterior
%    Variance VARF of the latent variable and returns the
%    posterior predictive mean EY and variance VARY of the
%    observations related to the latent variables. This subfunction
%    is needed when computing posterior predictive distributions for
%    future observations.
%

%
%  See also
%    GPLA_PRED, GPEP_PRED, GPMC_PRED
  
  N=sum(yt,2);
  S=10000;
  [ntest, nout]=size(yt);
  pi=zeros(ntest,nout);
  lpy=zeros(ntest,nout);
  Ey=zeros(ntest,nout);
  Vary=zeros(size(Varf));
  Ef=reshape(Ef(:),ntest,nout);
  [notused,notused,c] =size(Varf);
  if c>1
    mcmc=false;
  else
    mcmc=true;
    Varf=reshape(Varf(:), ntest, nout);
  end
  for i1=1:ntest
    if mcmc
      Sigm_tmp = (Varf(i1,:));
      f_star=bsxfun(@plus, Ef(i1,:), bsxfun(@times, sqrt(Sigm_tmp), ...
        randn(S,nout)));
    else
      Sigm_tmp=(Varf(:,:,i1)'+Varf(:,:,i1))./2;
      f_star=mvnrnd(Ef(i1,:), Sigm_tmp, S);
    end
    
    tmp = exp(f_star);
    tmp = tmp./(sum(tmp, 2)*ones(1,size(tmp,2)));
    
    if nargout > 1
        Ey(i1,:) = N(i1).*mean(tmp);
        for z1 = 1:nout;
          for z2 = 1:nout
            for z3=1:S
              Var_tmp(:,:,z3) = (diag(tmp(z3,:)) - tmp(z3,:)'*tmp(z3,:));
            end
            if mcmc
              Vary(i1+(0:nout-1)*ntest,:) = diag(N(i1).*mean(Var_tmp,3));
            else
              Vary(:,:,i1) = N(i1).*mean(Var_tmp,3);
            end
          end
        end
    end
    lpy=[];
    if ~isempty(yt)
      ytmp = repmat(yt(i1,:),S,1);
      lpy(i1,:) = log(mean( mnpdf(ytmp,tmp) ));
    end
  end
  lpy=lpy(:);
  Ey=Ey(:);
end

function p = lik_multinom_invlink(lik, f, z)
%LIK_MULTINOM_INVLINK Returns values of inverse link function
%             
%  Description 
%    P = LIK_MULTINOM_INVLINK(LIK, F) takes a likelihood structure LIK and
%    latent values F and returns the values of inverse link function P.
%    This subfunction is needed when using function gp_predprctmu.
%
%     See also
%     LIK_MULTINOM_LL, LIK_MULTINOM_PREDY
p = multinominv(f).*z;
end

function reclik = lik_multinom_recappend(reclik, ri, lik)
%RECAPPEND  Append the parameters to the record
%
%  Description 
%    RECLIK = LIK_MULTINOM_RECAPPEND(RECLIK, RI, LIK) takes a
%    likelihood record structure RECLIK, record index RI and
%    likelihood structure LIK with the current MCMC samples of
%    the parameters. Returns RECLIK which contains all the old
%    samples and the current samples from LIK. This subfunction 
%    is needed when using MCMC sampling (gp_mc).
% 
%  See also
%    GP_MC

  if nargin == 2
    reclik.type = 'Multinom';
    reclik.nondiagW = true;

    % Set the function handles
    reclik.fh.pak = @lik_multinom_pak;
    reclik.fh.unpak = @lik_multinom_unpak;
    reclik.fh.ll = @lik_multinom_ll;
    reclik.fh.llg = @lik_multinom_llg;    
    reclik.fh.llg2 = @lik_multinom_llg2;
    reclik.fh.llg3 = @lik_multinom_llg3;
    reclik.fh.predy = @lik_multinom_predy;
    reclik.fh.invlink = @lik_multinom_invlink;
    reclik.fh.recappend = @lik_multinom_recappend;
  end
  
end
