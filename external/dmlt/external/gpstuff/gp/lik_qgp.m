function lik = lik_qgp(varargin)
%LIK_QGP  Create a Quantile Gaussian Process likelihood (utility) structure
%
%  Description
%    LIK = LIK_QGP('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates a quantile gp likelihood structure in which the named
%    parameters have the specified values. Any unspecified
%    parameters are set to default values.
%
%    LIK = LIK_QGP(LIK,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a likelihood function structure with the named
%    parameters altered with the specified values.
%
%    Parameters for QGP likelihood function [default]
%      sigma2       - variance [0.1]
%      sigma2_prior - prior for sigma2 [prior_logunif]
%      quantile     - Quantile of interest [0.5]
%
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc. 
%
%    The likelihood is defined as follows:
%                            __ n
%      p(y|f, sigma2, tau) = || i=1 tau*(1-tau)/sigma*exp(-(y-f)/sigma*
%                                 (tau - I(t <= f)))
%    
%    where tau is the quantile of interest, sigma is the standard deviation
%    of the distribution and I(t <= f) = 1 if t <= f, 0 otherwise.
%
%    Note that because the form of the likelihood, second order derivatives
%    with respect to latent values are 0. Because this, EP should be used
%    instead of Laplace approximation.    
%
%  See also
%    GP_SET, PRIOR_*, LIK_*
%
%   References
%     Boukouvalas et al. (2012). Direct Gaussian Process Quantile Regression
%     Using Expectation Propagation. Appearing in Proceedings of the 29th
%     International Conference on Machine Learning, Edinburg, Scotland, UK,
%     2012.
%     
  
% Copyright (c) 2012 Ville Tolvanen

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_QGP';
  ip.addOptional('lik', [], @isstruct);
  ip.addParamValue('sigma2',0.1, @(x) isscalar(x) && x>0);
  ip.addParamValue('sigma2_prior',prior_logunif(), @(x) isstruct(x) || isempty(x));
  ip.addParamValue('quantile',0.5, @(x) isscalar(x) && x>0 && x<1);
  ip.parse(varargin{:});
  lik=ip.Results.lik;

  if isempty(lik)
    init=true;
    lik.type = 'QGP';
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'QGP')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end
  
  % Initialize parameters
  if init || ~ismember('sigma2',ip.UsingDefaults)
    lik.sigma2 = ip.Results.sigma2;
  end
  if init || ~ismember('quantile',ip.UsingDefaults)
    lik.quantile = ip.Results.quantile;
  end
  % Initialize prior structure
  if init
    lik.p=[];
  end
  if init || ~ismember('sigma2_prior',ip.UsingDefaults)
    lik.p.sigma2=ip.Results.sigma2_prior;
  end
  if init
    % Set the function handles to the subfunctions
    lik.fh.pak = @lik_qgp_pak;
    lik.fh.unpak = @lik_qgp_unpak;
    lik.fh.lp = @lik_qgp_lp;
    lik.fh.lpg = @lik_qgp_lpg;
    lik.fh.ll = @lik_qgp_ll;
    lik.fh.llg = @lik_qgp_llg;    
    lik.fh.llg2 = @lik_qgp_llg2;
    lik.fh.llg3 = @lik_qgp_llg3;
    lik.fh.tiltedMoments = @lik_qgp_tiltedMoments;
    lik.fh.siteDeriv = @lik_qgp_siteDeriv;
    lik.fh.predy = @lik_qgp_predy;
    lik.fh.invlink = @lik_qgp_invlink;
    lik.fh.recappend = @lik_qgp_recappend;
  end

end

function [w s] = lik_qgp_pak(lik)
%LIK_QGP_PAK  Combine likelihood parameters into one vector.
%
%  Description
%    W = LIK_QGP_PAK(LIK) takes a likelihood structure LIK
%    and combines the parameters into a single row vector W.
%    This is a mandatory subfunction used for example in 
%    energy and gradient computations.
%
%       w = [ log(lik.sigma2)
%             (hyperparameters of lik.magnSigma2)]'
%     
%  See also
%    LIK_QGP_UNPAK

  w = []; s = {};
  if ~isempty(lik.p.sigma2)
    w = [w log(lik.sigma2)];
    s = [s; 'log(qgp.sigma2)'];
    % Hyperparameters of sigma2
    [wh sh] = lik.p.sigma2.fh.pak(lik.p.sigma2);
    w = [w wh];
    s = [s; sh];
  end    

end

function [lik, w] = lik_qgp_unpak(lik, w)
%LIK_QGP_UNPAK  Extract likelihood parameters from the vector.
%
%  Description
%    W = LIK_QGP_UNPAK(W, LIK) takes a likelihood structure
%    LIK and extracts the parameters from the vector W to the LIK
%    structure. This is a mandatory subfunction used for example 
%    in energy and gradient computations.
%
%    Assignment is inverse of  
%       w = [ log(lik.sigma2)
%             (hyperparameters of lik.magnSigma2)]'
%
%  See also
%    LIK_QGP_PAK
  
  if ~isempty(lik.p.sigma2)
    lik.sigma2 = exp(w(1));
    w = w(2:end);
    
    % Hyperparameters of sigma2
    [p, w] = lik.p.sigma2.fh.unpak(lik.p.sigma2, w);
    lik.p.sigma2 = p;
  end
end

function lp = lik_qgp_lp(lik)
%LIK_QGP_LP  Evaluate the log prior of likelihood parameters
%
%  Description
%    LP = LIK_QGP_LP(LIK) takes a likelihood structure LIK and
%    returns log(p(th)), where th collects the parameters. This
%    subfunction is needed when there are likelihood parameters.
%
%  See also
%    LIK_QGP_PAK, LIK_QGP_UNPAK, LIK_QGP_G, GP_E

  lp = 0;

  if ~isempty(lik.p.sigma2)
    likp=lik.p;
    lp = likp.sigma2.fh.lp(lik.sigma2, likp.sigma2) + log(lik.sigma2);
  end
end

function lpg = lik_qgp_lpg(lik)
%LIK_QGP_LPG  Evaluate gradient of the log prior with respect
%                  to the parameters.
%
%  Description
%    LPG = LIK_QGP_LPG(LIK) takes a QGP likelihood
%    function structure LIK and returns LPG = d log (p(th))/dth,
%    where th is the vector of parameters. This subfunction is 
%    needed when there are likelihood parameters.
%
%  See also
%    LIK_QGP_PAK, LIK_QGP_UNPAK, LIK_QGP_E, GP_G

  lpg = [];

  if ~isempty(lik.p.sigma2)
    likp=lik.p;
    
    lpgs = likp.sigma2.fh.lpg(lik.sigma2, likp.sigma2);
    lpg = lpgs(1).*lik.sigma2 + 1;
    if length(lpgs) > 1
      lpg = [lpg lpgs(2:end)];
    end            
  end
end

function ll = lik_qgp_ll(lik, y, f, z)
%LIK_QGP_LL  Log likelihood
%
%  Description
%    LL = LIK_QGP_LL(LIK, Y, F, Z) takes a likelihood
%    structure LIK, observations Y and latent values F. 
%    Returns the log likelihood, log p(y|f,z). This subfunction 
%    is needed when using Laplace approximation or MCMC for 
%    inference with non-Gaussian likelihoods. This subfunction 
%    is also used in information criteria (DIC, WAIC) computations.
%
%  See also
%    LIK_QGP_LLG, LIK_QGP_LLG3, LIK_QGP_LLG2, GPLA_E
  
  tau=lik.quantile;
  sigma=sqrt(lik.sigma2);
  ll = sum(log(tau*(1-tau)/sigma) - (y-f)./sigma.*(tau-(y<=f)));
end

function llg = lik_qgp_llg(lik, y, f, param, z)
%LIK_QGP_LLG  Gradient of the log likelihood
%
%  Description 
%    LLG = LIK_QGP_LLG(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, observations Y and latent values F. Returns 
%    the gradient of the log likelihood with respect to PARAM. 
%    At the moment PARAM can be 'param' or 'latent'. This subfunction 
%    is needed when using Laplace approximation or MCMC for inference 
%    with non-Gaussian likelihoods.
%
%  See also
%    LIK_QGP_LL, LIK_QGP_LLG2, LIK_QGP_LLG3, GPLA_E

  
  tau=lik.quantile;
  sigma2=sqrt(lik.sigma2);
  switch param
    case 'param'      
      llg = sum(-1/(2.*sigma2) + (y-f)./(2.*sigma2^(3/2)).*(tau-(y<=f)));
      
      % correction for the log transformation
      llg = llg.*lik.sigma2;
    case 'latent'
      llg = (tau-(y<=f))/sqrt(sigma2);
  end
end

function llg2 = lik_qgp_llg2(lik, y, f, param, z)
%LIK_QGP_LLG2  Second gradients of the log likelihood
%
%  Description        
%    LLG2 = LIK_QGP_LLG2(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, observations Y and latent values F. Returns 
%    the Hessian of the log likelihood with respect to PARAM. 
%    At the moment PARAM can be 'param' or 'latent'. LLG2 is 
%    a vector with diagonal elements of the Hessian matrix 
%    (off diagonals are zero). This subfunction is needed 
%    when using Laplace approximation or EP for inference 
%    with non-Gaussian likelihoods.
%
%  See also
%    LIK_QGP_LL, LIK_QGP_LLG, LIK_QGP_LLG3, GPLA_E

  
  tau=lik.quantile;
  sigma2=lik.sigma2;
  switch param
    case 'param'
      llg2 = sum(1/(2*sigma2^2) - 3.*(tau-(y<=f)).*(y-f)./(4.*sigma2^(5/2)));
      
      % correction due to the log transformation
      llg2 = llg2.*lik.sigma2;
    case 'latent'
      llg2 = zeros(size(f));
    case 'latent+param'
      llg2 = -(tau-(y<=f))./(2*sigma2^(3/2));
      
      % correction due to the log transformation
      llg2 = llg2.*lik.disper;
  end
end    

function llg3 = lik_qgp_llg3(lik, y, f, param, z)
%LIK_QGP_LLG3  Third gradients of the log likelihood
%
%  Description
%    LLG3 = LIK_QGP_LLG3(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, observations Y and latent values F and 
%    returns the third gradients of the log likelihood with 
%    respect to PARAM. At the moment PARAM can be 'param' or 
%    'latent'. LLG3 is a vector with third gradients. This 
%    subfunction is needed when using Laplace approximation for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_QGP_LL, LIK_QGP_LLG, LIK_QGP_LLG2, GPLA_E, GPLA_G

  tau=lik.quantile;
  sigma2=lik.sigma2;
  switch param
    case 'param'
      llg3 = sum(-1/sigma2^3 + 15.*(tau-(y<=f)).*(y-f)./(8.*sigma2^(7/2)));
    case 'latent'
      llg3 = 0;
    case 'latent2+param'
      llg3 = 0;
      
      % correction due to the log transformation
      llg3 = llg3.*lik.sigma2;
  end
end

function [logM_0, m_1, sigm2hati1] = lik_qgp_tiltedMoments(lik, y, i1, sigm2_i, myy_i, z)
%LIK_QGP_TILTEDMOMENTS  Returns the marginal moments for EP algorithm
%
%  Description
%    [M_0, M_1, M2] = LIK_QGP_TILTEDMOMENTS(LIK, Y, I, S2,
%    MYY, Z) takes a likelihood structure LIK, observations
%    Y, index I and cavity variance S2 and mean MYY. Returns 
%    the zeroth moment M_0, mean M_1 and variance M_2 of the 
%    posterior marginal (see Rasmussen and Williams (2006): 
%    Gaussian processes for Machine Learning, page 55). This 
%    subfunction is needed when using EP for inference with 
%    non-Gaussian likelihoods.
%
%  See also
%    GPEP_E
  
  yy = y(i1);
  sigma2 = lik.sigma2;
  tau=lik.quantile;
  logM_0=zeros(size(yy));
  m_1=zeros(size(yy));
  sigm2hati1=zeros(size(yy));
  
  for i=1:length(i1)
    % get a function handle of an unnormalized tilted distribution
    % (likelihood * cavity = Quantile-GP * Gaussian)
    % and useful integration limits
    [tf,minf,maxf]=init_qgp_norm(yy(i),myy_i(i),sigm2_i(i),sigma2,tau);
    
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
        error('lik_qgp_tilted_moments: sigm2hati1 >= sigm2_i');
      end
    end
    logM_0(i) = log(m_0);
  end
end

function [g_i] = lik_qgp_siteDeriv(lik, y, i1, sigm2_i, myy_i, z)
%LIK_QGP_SITEDERIV  Evaluate the expectation of the gradient
%                      of the log likelihood term with respect
%                      to the likelihood parameters for EP 
%
%  Description [M_0, M_1, M2] =
%    LIK_QGP_SITEDERIV(LIK, Y, I, S2, MYY, Z) takes a
%    likelihood structure LIK, observations Y, index I 
%    and cavity variance S2 and mean MYY. Returns E_f 
%    [d log p(y_i|f_i) /d a], where a is the likelihood 
%    parameter and the expectation is over the marginal posterior.
%    This term is needed when evaluating the gradients of 
%    the marginal likelihood estimate Z_EP with respect to 
%    the likelihood parameters (see Seeger (2008):
%    Expectation propagation for exponential families).This 
%    subfunction is needed when using EP for inference with 
%    non-Gaussian likelihoods and there are likelihood parameters.
%
%  See also
%    GPEP_G


  yy = y(i1);
  sigma2=lik.sigma2;
  tau=lik.quantile;
  
  % get a function handle of an unnormalized tilted distribution 
  % (likelihood * cavity = Quantile-GP * Gaussian)
  % and useful integration limits
  [tf,minf,maxf]=init_qgp_norm(yy,myy_i,sigm2_i,sigma2,tau);
  % additionally get function handle for the derivative
  td = @deriv;
  
  % Integrate with quadgk
  [m_0, fhncnt] = quadgk(tf, minf, maxf);
  [g_i, fhncnt] = quadgk(@(f) td(f).*tf(f)./m_0, minf, maxf);
  g_i = g_i.*sigma2;

  function g = deriv(f)

    g = -1/(2.*sigma2) + (yy-f)./(2.*sigma2^(3/2)).*(tau-(yy<=f));
    
  end
end

function [lpy, Ey, Vary] = lik_qgp_predy(lik, Ef, Varf, yt, zt)
%LIK_QGP_PREDY  Returns the predictive mean, variance and density of y
%
%  Description  
%    LPY = LIK_QGP_PREDY(LIK, EF, VARF YT, ZT)
%    Returns logarithm of the predictive density PY of YT, that is 
%        p(yt | zt) = \int p(yt | f, zt) p(f|y) df.
%    This subfunction is needed when computing posterior predictive 
%    distributions for future observations.
%
%    [LPY, EY, VARY] = LIK_QGP_PREDY(LIK, EF, VARF) takes a
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


  sigma2=lik.sigma2;
  tau=lik.quantile;
  
  Ey=[];
  Vary=[];
  
  % Evaluate the posterior predictive densities of the given observations
  lpy = zeros(length(yt),1);
  for i1=1:length(yt)
    % get a function handle of the likelihood times posterior
    % (likelihood * posterior = Quantile-GP * Gaussian)
    % and useful integration limits
    [pdf,minf,maxf]=init_qgp_norm(...
      yt(i1),Ef(i1),Varf(i1),sigma2, tau);
    % integrate over the f to get posterior predictive distribution
    lpy(i1) = log(quadgk(pdf, minf, maxf));
  end
end


function [df,minf,maxf] = init_qgp_norm(yy,myy_i,sigm2_i,sigma2,tau)
%INIT_QGP_NORM
%
%  Description
%    Return function handle to a function evaluating
%    Quantile-GP * Gaussian which is used for evaluating
%    (likelihood * cavity) or (likelihood * posterior) Return
%    also useful limits for integration. This is private function
%    for lik_qgp. This subfunction is needed by subfunctions
%    tiltedMoments, siteDeriv and predy.
%  
%  See also
%    LIK_QGP_TILTEDMOMENTS, LIK_QGP_SITEDERIV,
%    LIK_QGP_PREDY
  
  sigma=sqrt(sigma2);
% avoid repetitive evaluation of constant part
  ldconst = log(tau*(1-tau)/sigma) ...
            - log(sigm2_i)/2 - log(2*pi)/2;
  % Create function handle for the function to be integrated
  df = @qgp_norm;
  % use log to avoid underflow, and derivates for faster search
  ld = @log_qgp_norm;
  ldg = @log_qgp_norm_g;
%   ldg2 = @log_qgp_norm_g2;

  % Set the limits for integration
  % Quantile-GP likelihood is log-concave so the qgp_norm
  % function is unimodal, which makes things easier
  if yy==0
    % with yy==0, the mode of the likelihood is not defined
    % use the mode of the Gaussian (cavity or posterior) as a first guess
    modef = myy_i;
  else
    % use precision weighted mean of the Gaussian approximation
    % of the Quantile-GP likelihood and Gaussian
    modef = (myy_i/sigm2_i + yy/sigma2)/(1/sigm2_i + 1/sigma2);
  end
  % find the mode of the integrand using Newton iterations
  % few iterations is enough, since the first guess in the right direction
  niter=8;       % number of Newton iterations 
  
  minf=modef-6*sigm2_i;
  while ldg(minf) < 0
    minf=minf-2*sigm2_i;
  end
  maxf=modef+6*sigm2_i;
  while ldg(maxf) > 0
    maxf=maxf+2*sigm2_i;
  end
  for ni=1:niter
%     h=ldg2(modef);
    modef=0.5*(minf+maxf);
    if ldg(modef) < 0
      maxf=modef;
    else
      minf=modef;
    end
  end
  % integrand limits based on Gaussian approximation at mode
  minf=modef-6*sqrt(sigm2_i);
  maxf=modef+6*sqrt(sigm2_i);
  modeld=ld(modef);
  iter=0;
  % check that density at end points is low enough
  lddiff=20; % min difference in log-density between mode and end-points
  minld=ld(minf);
  step=1;
  while minld>(modeld-lddiff)
    minf=minf-step*sqrt(sigm2_i);
    minld=ld(minf);
    iter=iter+1;
    step=step*2;
    if iter>100
      error(['lik_qgp -> init_qgp_norm: ' ...
             'integration interval minimun not found ' ...
             'even after looking hard!'])
    end
  end
  maxld=ld(maxf);
  step=1;
  while maxld>(modeld-lddiff)
    maxf=maxf+step*sqrt(sigm2_i);
    maxld=ld(maxf);
    iter=iter+1;
    step=step*2;
    if iter>100
      error(['lik_qgp -> init_qgp_norm: ' ...
             'integration interval maximun not found ' ...
             'even after looking hard!'])
    end
  end
  
  function integrand = qgp_norm(f)
  % Quantile-GP * Gaussian
    integrand = exp(ldconst ...
                    -(yy-f)./sqrt(sigma2).*(tau-(yy<=f)) ...
                    -0.5*(f-myy_i).^2./sigm2_i);
  end
  
  function log_int = log_qgp_norm(f)
  % log(Quantile-GP * Gaussian)
  % log_qgp_norm is used to avoid underflow when searching
  % integration interval
    log_int = ldconst...
              -(yy-f)./sqrt(sigma2).*(tau-(yy<=f)) ...
              -0.5*(f-myy_i).^2./sigm2_i;
  end
  
  function g = log_qgp_norm_g(f)
  % d/df log(Quantile-GP * Gaussian)
  % derivative of log_qgp_norm
    g = (tau-(yy<=f))/sqrt(sigma2) ...
        + (myy_i - f)./sigm2_i;
  end
  
  
end

function mu = lik_qgp_invlink(lik, f, z)
%LIK_QGP_INVLINK  Returns values of inverse link function
%             
%  Description 
%    MU = LIK_QGP_INVLINK(LIK, F) takes a likelihood structure LIK and
%    latent values F and returns the values MU of inverse link function.
%    This subfunction is needed when using function gp_predprctmu.
%
%     See also
%     LIK_QGP_LL, LIK_QGP_PREDY
  
  mu = f;
end

function reclik = lik_qgp_recappend(reclik, ri, lik)
%RECAPPEND  Append the parameters to the record
%
%  Description 
%    RECLIK = LIK_QGP_RECAPPEND(RECLIK, RI, LIK) takes a
%    likelihood record structure RECLIK, record index RI and
%    likelihood structure LIK with the current MCMC samples of
%    the parameters. Returns RECLIK which contains all the old
%    samples and the current samples from LIK.  This subfunction
%    is needed when using MCMC sampling (gp_mc).
% 
%  See also
%    GP_MC

  if nargin == 2
    % Initialize the record
    reclik.type = 'Quantile-GP';

    % Initialize parameter
    reclik.sigma2 = [];

    % Set the function handles
    reclik.fh.pak = @lik_qgp_pak;
    reclik.fh.unpak = @lik_qgp_unpak;
    reclik.fh.lp = @lik_qgp_lp;
    reclik.fh.lpg = @lik_qgp_lpg;
    reclik.fh.ll = @lik_qgp_ll;
    reclik.fh.llg = @lik_qgp_llg;    
    reclik.fh.llg2 = @lik_qgp_llg2;
    reclik.fh.llg3 = @lik_qgp_llg3;
    reclik.fh.tiltedMoments = @lik_qgp_tiltedMoments;
    reclik.fh.predy = @lik_qgp_predy;
    reclik.fh.invlink = @lik_qgp_invlink;
    reclik.fh.recappend = @lik_qgp_recappend;
    reclik.p=[];
    reclik.p.sigma2=[];
    if ~isempty(ri.p.sigma2)
      reclik.p.sigma2 = ri.p.sigma2;
    end
  else
    
    % Append to the record
    reclik.sigma2(ri,:)=lik.sigma2;
    if ~isempty(lik.p)
      reclik.p.sigma2 = lik.p.sigma2.fh.recappend(reclik.p.sigma2, ri, lik.p.sigma2);
    end
  end
end
