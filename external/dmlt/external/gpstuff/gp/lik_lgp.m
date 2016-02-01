function lik = lik_lgp(varargin)
%LIK_LGP  Create a logistic Gaussian process likelihood structure 
%
%  Description
%    LIK = LIK_LGP creates a logistic Gaussian process likelihood structure
%
%    The likelihood is defined as follows:
%               __ n
%      p(y|f) = || i=1 exp(f_i) / Sum_{j=1}^n exp(f_j),
%
%      where f contains latent values.
%
%  See also
%    LGPDENS, GP_SET, LIK_*
%
  
% Copyright (c) 2011 Jaakko Riihimäki and Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_LGP';
  ip.addOptional('lik', [], @isstruct);
  ip.parse(varargin{:});
  lik=ip.Results.lik;

  if isempty(lik)
    init=true;
    lik.type = 'LGP';
    lik.nondiagW = true;
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'LGP')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end

  if init
    % Set the function handles to the subfunctions
    lik.fh.pak = @lik_lgp_pak;
    lik.fh.unpak = @lik_lgp_unpak;
    lik.fh.ll = @lik_lgp_ll;
    lik.fh.llg = @lik_lgp_llg;    
    lik.fh.llg2 = @lik_lgp_llg2;
    lik.fh.llg3 = @lik_lgp_llg3;
    lik.fh.tiltedMoments = @lik_lgp_tiltedMoments;
    lik.fh.predy = @lik_lgp_predy;
    lik.fh.invlink = @lik_lgp_invlink;
    lik.fh.recappend = @lik_lgp_recappend;
  end

end

function [w,s] = lik_lgp_pak(lik)
%LIK_LGP_PAK  Combine likelihood parameters into one vector.
%
%  Description 
%    W = LIK_LGP_PAK(LIK) takes a likelihood structure LIK
%    and returns an empty verctor W. If LGP likelihood had
%    parameters this would combine them into a single row vector
%    W (see e.g. lik_negbin). This is a mandatory subfunction 
%    used for example in energy and gradient computations.
%     
%  See also
%    LIK_LGP_UNPAK, GP_PAK

  w = []; s = {};
end


function [lik, w] = lik_lgp_unpak(lik, w)
%LIK_LGP_UNPAK  Extract likelihood parameters from the vector.
%
%  Description
%    W = LIK_LGP_UNPAK(W, LIK) Doesn't do anything.
%
%    If LGP likelihood had parameters this would extract them
%    parameters from the vector W to the LIK structure. This 
%    is a mandatory subfunction used for example in energy 
%    and gradient computations.
%     
%
%  See also
%    LIK_LGP_PAK, GP_UNPAK

  lik=lik;
  w=w;
  
end


function logLik = lik_lgp_ll(lik, y, f, z)
%LIK_LGP_LL    Log likelihood
%
%  Description
%    E = LIK_LGP_LL(LIK, Y, F, Z) takes a likelihood data
%    structure LIK, incedence counts Y, expected counts Z, and
%    latent values F. Returns the log likelihood, log p(y|f,z).
%    This subfunction is needed when using Laplace approximation 
%    or MCMC for inference with non-Gaussian likelihoods. This 
%    subfunction is also used in information criteria (DIC, WAIC) 
%    computations.
%
%  See also
%    LIK_LGP_LLG, LIK_LGP_LLG3, LIK_LGP_LLG2, GPLA_E

  n=sum(y);
  qj=exp(f);
  logLik = sum(f.*y)-n*log(sum(qj));
end


function deriv = lik_lgp_llg(lik, y, f, param, z)
%LIK_LGP_LLG    Gradient of the log likelihood
%
%  Description 
%    G = LIK_LGP_LLG(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, incedence counts Y, expected counts Z
%    and latent values F. Returns the gradient of the log
%    likelihood with respect to PARAM. At the moment PARAM can be
%    'param' or 'latent'. This subfunction is needed when using Laplace 
%    approximation or MCMC for inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_LGP_LL, LIK_LGP_LLG2, LIK_LGP_LLG3, GPLA_E
  
  switch param
    case 'latent'
      n=sum(y);
      qj=exp(f);
      pj=qj./sum(qj);
      deriv=y-n*pj;
  end
end


function g2 = lik_lgp_llg2(lik, y, f, param, z)
%function g2 = lik_lgp_llg2(lik, y, f, param, z)
%LIK_LGP_LLG2  Second gradients of the log likelihood
%
%  Description        
%    G2 = LIK_LGP_LLG2(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, incedence counts Y, expected counts Z,
%    and latent values F. Returns the Hessian of the log
%    likelihood with respect to PARAM. At the moment PARAM can be
%    only 'latent'. G2 is a vector with diagonal elements of the
%    Hessian matrix (off diagonals are zero). This subfunction 
%    is needed when using Laplace approximation or EP for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_LGP_LL, LIK_LGP_LLG, LIK_LGP_LLG3, GPLA_E

  switch param
    case 'latent'
      qj=exp(f);
      
      % g2 is not the second gradient of the log likelihood but only a
      % vector to form the exact gradient term in gpla_nd_e, gpla_nd_g and
      % gpla_nd_pred functions
      g2=qj./sum(qj);
  end
end    

function g3 = lik_lgp_llg3(lik, y, f, param, z)
%LIK_LGP_LLG3  Third gradients of the log likelihood
%
%  Description
%    G3 = LIK_LGP_LLG3(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, incedence counts Y, expected counts Z
%    and latent values F and returns the third gradients of the
%    log likelihood with respect to PARAM. At the moment PARAM
%    can be only 'latent'. G3 is a vector with third gradients.
%    This subfunction is needed when using Laplace approximation 
%    for inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_LGP_LL, LIK_LGP_LLG, LIK_LGP_LLG2, GPLA_E, GPLA_G
  
  switch param
    case 'latent'
      qj=exp(f);
      
      % g3 is not the third gradient of the log likelihood but only a
      % vector to form the exact gradient term in gpla_nd_e, gpla_nd_g and
      % gpla_nd_pred functions
      g3=qj./sum(qj);
      
      %n=sum(y);
      %nf=size(f,1);
      %g3d=zeros(nf,nf);
      %for i1=1:nf
      %  g3dtmp=-g3*g3(i1);
      %  g3dtmp(i1)=g3dtmp(i1)+g3(i1);
      %  g3d(:,i1)=g3dtmp;
      %  %g3i1= n*(-diag(g3d(:,i1)) + bsxfun(@times,g3,g3d(:,i1)') + bsxfun(@times,g3d(:,i1),g3'));
      %end
  end
end

function [logM_0, m_1, sigm2hati1] = lik_lgp_tiltedMoments(lik, y, i1, sigm2_i, myy_i, z)
%LIK_LGP_TILTEDMOMENTS  Returns the marginal moments for EP algorithm
%
%  Description
%    [M_0, M_1, M2] = LIK_LGP_TILTEDMOMENTS(LIK, Y, I, S2,
%    MYY, Z) takes a likelihood structure LIK, incedence counts
%    Y, expected counts Z, index I and cavity variance S2 and
%    mean MYY. Returns the zeroth moment M_0, mean M_1 and
%    variance M_2 of the posterior marginal (see Rasmussen and
%    Williams (2006): Gaussian processes for Machine Learning,
%    page 55). This subfunction is needed when using EP for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    GPEP_E

  
  if isempty(z)
    error(['lik_lgp -> lik_lgp_tiltedMoments: missing z!'... 
           'LGP likelihood needs the expected number of             '...
           'occurrences as an extra input z. See, for                   '...
           'example, lik_lgp and gpla_e.                        ']);
  end
  
  yy = y(i1);
  avgE = z(i1);
  logM_0=zeros(size(yy));
  m_1=zeros(size(yy));
  sigm2hati1=zeros(size(yy));  
  
  for i=1:length(i1)
    % get a function handle of an unnormalized tilted distribution
    % (likelihood * cavity = Negative-binomial * Gaussian)
    % and useful integration limits
    [tf,minf,maxf]=init_lgp_norm(yy(i),myy_i(i),sigm2_i(i),avgE(i));
    
    % Integrate with quadrature
    RTOL = 1.e-6;
    ATOL = 1.e-10;
    [m_0, m_1(i), m_2] = quad_moments(tf, minf, maxf, RTOL, ATOL);
    sigm2hati1(i) = m_2 - m_1(i).^2;
    
    % If the second central moment is less than cavity variance
    % integrate more precisely. Theoretically for log-concave
    % likelihood should be sigm2hati1 < sigm2_i.
    if sigm2hati1(i) >= sigm2_i(i)
      ATOL = ATOL.^2;
      RTOL = RTOL.^2;
      [m_0, m_1(i), m_2] = quad_moments(tf, minf, maxf, RTOL, ATOL);
      sigm2hati1(i) = m_2 - m_1(i).^2;
      if sigm2hati1(i) >= sigm2_i(i)
        error('lik_lgp_tilted_moments: sigm2hati1 >= sigm2_i');
      end
    end
    logM_0(i) = log(m_0);
  end
end


function [lpy, Ey, Vary] = lik_lgp_predy(lik, Ef, Varf, yt, zt)
%LIK_LGP_PREDY    Returns the predictive mean, variance and density of y
%
%  Description  
%    LPY = LIK_LGP_PREDY(LIK, EF, VARF YT, ZT)
%    Returns also the predictive density of YT, that is 
%        p(yt | y,zt) = \int p(yt | f, zt) p(f|y) df.
%    This requires also the incedence counts YT, expected counts ZT.
%    This subfunction is needed when computing posterior predictive 
%    distributions for future observations.
%
%    [LPY, EY, VARY] = LIK_LGP_PREDY(LIK, EF, VARF) takes a
%    likelihood structure LIK, posterior mean EF and posterior
%    Variance VARF of the latent variable and returns the
%    posterior predictive mean EY and variance VARY of the
%    observations related to the latent variables. This 
%    subfunction is needed when computing posterior predictive 
%    distributions for future observations.
%        

%
%  See also 
%    GPLA_PRED, GPEP_PRED, GPMC_PRED

  if isempty(zt)
    error(['lik_lgp -> lik_lgp_predy: missing zt!'... 
           'LGP likelihood needs the expected number of     '...
           'occurrences as an extra input zt. See, for           '...
           'example, lik_lgp and gpla_e.                ']);
  end
  
  avgE = zt;
  lpy = zeros(size(Ef));
  Ey = zeros(size(Ef));
  EVary = zeros(size(Ef));
  VarEy = zeros(size(Ef)); 
  
  if nargout > 1
      % Evaluate Ey and Vary
      for i1=1:length(Ef)
        %%% With quadrature
        myy_i = Ef(i1);
        sigm_i = sqrt(Varf(i1));
        minf=myy_i-6*sigm_i;
        maxf=myy_i+6*sigm_i;

        F = @(f) exp(log(avgE(i1))+f+norm_lpdf(f,myy_i,sigm_i));
        Ey(i1) = quadgk(F,minf,maxf);

        EVary(i1) = Ey(i1);

        F3 = @(f) exp(2*log(avgE(i1))+2*f+norm_lpdf(f,myy_i,sigm_i));
        VarEy(i1) = quadgk(F3,minf,maxf) - Ey(i1).^2;
      end
      Vary = EVary + VarEy;
  end

  % Evaluate the posterior predictive densities of the given observations
  for i1=1:length(Ef)
    % get a function handle of the likelihood times posterior
    % (likelihood * posterior = LGP * Gaussian)
    % and useful integration limits
    [pdf,minf,maxf]=init_lgp_norm(...
      yt(i1),Ef(i1),Varf(i1),avgE(i1));
    % integrate over the f to get posterior predictive distribution
    lpy(i1) = log(quadgk(pdf, minf, maxf));
  end
end

function [df,minf,maxf] = init_lgp_norm(yy,myy_i,sigm2_i,avgE)
%INIT_LGP_NORM
%
%  Description
%    Return function handle to a function evaluating LGP *
%    Gaussian which is used for evaluating (likelihood * cavity)
%    or (likelihood * posterior) Return also useful limits for
%    integration. This is private function for lik_lgp. This 
%    subfunction is needed by sufunctions tiltedMoments, siteDeriv 
%    and predy.
%  
%  See also
%    LIK_LGP_TILTEDMOMENTS, LIK_LGP_PREDY
  
% avoid repetitive evaluation of constant part
  ldconst = -gammaln(yy+1) - log(sigm2_i)/2 - log(2*pi)/2;
  
  % Create function handle for the function to be integrated
  df = @lgp_norm;
  % use log to avoid underflow, and derivates for faster search
  ld = @log_lgp_norm;
  ldg = @log_lgp_norm_g;
  ldg2 = @log_lgp_norm_g2;

  % Set the limits for integration
  % LGP likelihood is log-concave so the lgp_norm
  % function is unimodal, which makes things easier
  if yy==0
    % with yy==0, the mode of the likelihood is not defined
    % use the mode of the Gaussian (cavity or posterior) as a first guess
    modef = myy_i;
  else
    % use precision weighted mean of the Gaussian approximation
    % of the LGP likelihood and Gaussian
    mu=log(yy/avgE);
    s2=1./(yy+1./sigm2_i);
    modef = (myy_i/sigm2_i + mu/s2)/(1/sigm2_i + 1/s2);
  end
  % find the mode of the integrand using Newton iterations
  % few iterations is enough, since the first guess in the right direction
  niter=3;       % number of Newton iterations
  mindelta=1e-6; % tolerance in stopping Newton iterations
  for ni=1:niter
    g=ldg(modef);
    h=ldg2(modef);
    delta=-g/h;
    modef=modef+delta;
    if abs(delta)<mindelta
      break
    end
  end
  % integrand limits based on Gaussian approximation at mode
  modes=sqrt(-1/h);
  minf=modef-8*modes;
  maxf=modef+8*modes;
  modeld=ld(modef);
  iter=0;
  % check that density at end points is low enough
  lddiff=20; % min difference in log-density between mode and end-points
  minld=ld(minf);
  step=1;
  while minld>(modeld-lddiff)
    minf=minf-step*modes;
    minld=ld(minf);
    iter=iter+1;
    step=step*2;
    if iter>100
      error(['lik_lgp -> init_lgp_norm: ' ...
             'integration interval minimun not found ' ...
             'even after looking hard!'])
    end
  end
  maxld=ld(maxf);
  step=1;
  while maxld>(modeld-lddiff)
    maxf=maxf+step*modes;
    maxld=ld(maxf);
    iter=iter+1;
    step=step*2;
    if iter>100
      error(['lik_lgp -> init_lgp_norm: ' ...
             'integration interval maximun not found ' ...
             'even after looking hard!'])
    end
  end

  function integrand = lgp_norm(f)
  % LGP * Gaussian
    mu = avgE.*exp(f);
    integrand = exp(ldconst ...
                    -mu+yy.*log(mu) ...
                    -0.5*(f-myy_i).^2./sigm2_i);
  end

  function log_int = log_lgp_norm(f)
  % log(LGP * Gaussian)
  % log_lgp_norm is used to avoid underflow when searching
  % integration interval
    mu = avgE.*exp(f);
    log_int = ldconst ...
              -mu+yy.*log(mu) ...
              -0.5*(f-myy_i).^2./sigm2_i;
  end

  function g = log_lgp_norm_g(f)
  % d/df log(LGP * Gaussian)
  % derivative of log_lgp_norm
    mu = avgE.*exp(f);
    g = -mu+yy...
        + (myy_i - f)./sigm2_i;
  end

  function g2 = log_lgp_norm_g2(f)
  % d^2/df^2 log(LGP * Gaussian)
  % second derivate of log_lgp_norm
    mu = avgE.*exp(f);
    g2 = -mu...
         -1/sigm2_i;
  end

end

function mu = lik_lgp_invlink(lik, f, z)
%LIK_LGP_INVLINK  Returns values of inverse link function
%             
%  Description 
%    P = LIK_LGP_INVLINK(LIK, F) takes a likelihood structure LIK and
%    latent values F and returns the values MU of inverse link function.
%    This subfunction is needed when using function gp_predprctmu.
%
%     See also
%     LIK_LGP_LL, LIK_LGP_PREDY
  
  mu = z.*exp(f);
end

function reclik = lik_lgp_recappend(reclik, ri, lik)
%RECAPPEND  Append the parameters to the record
%
%  Description 
%    RECLIK = LIK_LGP_RECAPPEND(RECLIK, RI, LIK) takes a
%    likelihood record structure RECLIK, record index RI and
%    likelihood structure LIK with the current MCMC samples of
%    the parameters. Returns RECLIK which contains all the old
%    samples and the current samples from LIK. This subfunction 
%    is needed when using MCMC sampling (gp_mc).
% 
%  See also
%    GP_MC

  if nargin == 2
    reclik.type = 'LGP';

    % Set the function handles
    reclik.fh.pak = @lik_lgp_pak;
    reclik.fh.unpak = @lik_lgp_unpak;
    reclik.fh.ll = @lik_lgp_ll;
    reclik.fh.llg = @lik_lgp_llg;    
    reclik.fh.llg2 = @lik_lgp_llg2;
    reclik.fh.llg3 = @lik_lgp_llg3;
    reclik.fh.tiltedMoments = @lik_lgp_tiltedMoments;
    reclik.fh.predy = @lik_lgp_predy;
    reclik.fh.invlink = @lik_lgp_invlink;
    reclik.fh.recappend = @lik_lgp_recappend;
    return
  end
end
