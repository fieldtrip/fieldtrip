function lik = lik_weibull(varargin)
%LIK_WEIBULL    Create a right censored Weibull likelihood structure 
%
%  Description
%    LIK = LIK_WEIBULL('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates a likelihood structure for right censored Weibull
%    survival model in which the named parameters have the
%    specified values. Any unspecified parameters are set to
%    default values.
%  
%    LIK = LIK_WEIBULL(LIK,'PARAM1',VALUE1,'PARAM2,VALUE2,...)
%    modify a likelihood structure with the named parameters
%    altered with the specified values.
%
%    Parameters for Weibull likelihood [default]
%      shape       - shape parameter r [1]
%      shape_prior - prior for shape [prior_logunif]
%  
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc.
%
%    The likelihood is defined as follows:
%                  __ n
%      p(y|f, z) = || i=1 [ r^(1-z_i) exp( (1-z_i)*(-f_i)
%                           +(1-z_i)*(r-1)*log(y_i)
%                           -exp(-f_i)*y_i^r) ]
%
%    where r is the shape parameter of Weibull distribution.
%    z is a vector of censoring indicators with z = 0 for uncensored event
%    and z = 1 for right censored event. 
%
%    When using the Weibull likelihood you need to give the vector z
%    as an extra parameter to each function that requires also y. 
%    For example, you should call gpla_e as follows: gpla_e(w, gp,
%    x, y, 'z', z)
%
%  See also
%    GP_SET, LIK_*, PRIOR_*
%

% Copyright (c) 2011 Jaakko Riihimäki
% Copyright (c) 2011 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_WEIBULL';
  ip.addOptional('lik', [], @isstruct);
  ip.addParamValue('shape',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('shape_prior',prior_logunif(), @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  lik=ip.Results.lik;
  
  if isempty(lik)
    init=true;
    lik.type = 'Weibull';
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'Weibull')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('shape',ip.UsingDefaults)
    lik.shape = ip.Results.shape;
  end
  % Initialize prior structure
  if init
    lik.p=[];
  end
  if init || ~ismember('shape_prior',ip.UsingDefaults)
    lik.p.shape=ip.Results.shape_prior;
  end
  
  if init
    % Set the function handles to the subfunctions
    lik.fh.pak = @lik_weibull_pak;
    lik.fh.unpak = @lik_weibull_unpak;
    lik.fh.lp = @lik_weibull_lp;
    lik.fh.lpg = @lik_weibull_lpg;
    lik.fh.ll = @lik_weibull_ll;
    lik.fh.llg = @lik_weibull_llg;    
    lik.fh.llg2 = @lik_weibull_llg2;
    lik.fh.llg3 = @lik_weibull_llg3;
    lik.fh.tiltedMoments = @lik_weibull_tiltedMoments;
    lik.fh.siteDeriv = @lik_weibull_siteDeriv;
    lik.fh.invlink = @lik_weibull_invlink;
    lik.fh.predy = @lik_weibull_predy;
    lik.fh.predcdf = @lik_weibull_predcdf;
    lik.fh.recappend = @lik_weibull_recappend;
  end

end

function [w,s] = lik_weibull_pak(lik)
%LIK_WEIBULL_PAK  Combine likelihood parameters into one vector.
%
%  Description 
%    W = LIK_WEIBULL_PAK(LIK) takes a likelihood structure LIK and
%    combines the parameters into a single row vector W. This is a 
%    mandatory subfunction used for example in energy and gradient 
%    computations.
%     
%       w = log(lik.shape)
%
%   See also
%   LIK_WEIBULL_UNPAK, GP_PAK
  
  w=[];s={};
  if ~isempty(lik.p.shape)
    w = log(lik.shape);
    s = [s; 'log(weibull.shape)'];
    [wh sh] = lik.p.shape.fh.pak(lik.p.shape);
    w = [w wh];
    s = [s; sh];
  end
end


function [lik, w] = lik_weibull_unpak(lik, w)
%LIK_WEIBULL_UNPAK  Extract likelihood parameters from the vector.
%
%  Description
%    [LIK, W] = LIK_WEIBULL_UNPAK(W, LIK) takes a likelihood
%    structure LIK and extracts the parameters from the vector W
%    to the LIK structure. This is a mandatory subfunction used 
%    for example in energy and gradient computations.
%     
%   Assignment is inverse of  
%       w = log(lik.shape)
%
%   See also
%   LIK_WEIBULL_PAK, GP_UNPAK

  if ~isempty(lik.p.shape)
    lik.shape = exp(w(1));
    w = w(2:end);
    [p, w] = lik.p.shape.fh.unpak(lik.p.shape, w);
    lik.p.shape = p;
  end
end


function lp = lik_weibull_lp(lik, varargin)
%LIK_WEIBULL_LP  log(prior) of the likelihood parameters
%
%  Description
%    LP = LIK_WEIBULL_LP(LIK) takes a likelihood structure LIK and
%    returns log(p(th)), where th collects the parameters. This 
%    subfunction is needed when there are likelihood parameters.
%
%  See also
%    LIK_WEIBULL_LLG, LIK_WEIBULL_LLG3, LIK_WEIBULL_LLG2, GPLA_E
  

% If prior for shape parameter, add its contribution
  lp=0;
  if ~isempty(lik.p.shape)
    lp = lik.p.shape.fh.lp(lik.shape, lik.p.shape) +log(lik.shape);
  end
  
end


function lpg = lik_weibull_lpg(lik)
%LIK_WEIBULL_LPG  d log(prior)/dth of the likelihood 
%                parameters th
%
%  Description
%    E = LIK_WEIBULL_LPG(LIK) takes a likelihood structure LIK and
%    returns d log(p(th))/dth, where th collects the parameters.
%    This subfunction is needed when there are likelihood parameters.
%
%  See also
%    LIK_WEIBULL_LLG, LIK_WEIBULL_LLG3, LIK_WEIBULL_LLG2, GPLA_G
  
  lpg=[];
  if ~isempty(lik.p.shape)            
    % Evaluate the gprior with respect to shape
    ggs = lik.p.shape.fh.lpg(lik.shape, lik.p.shape);
    lpg = ggs(1).*lik.shape + 1;
    if length(ggs) > 1
      lpg = [lpg ggs(2:end)];
    end
  end
end  

function ll = lik_weibull_ll(lik, y, f, z)
%LIK_WEIBULL_LL  Log likelihood
%
%  Description
%    LL = LIK_WEIBULL_LL(LIK, Y, F, Z) takes a likelihood
%    structure LIK, survival times Y, censoring indicators Z, and
%    latent values F. Returns the log likelihood, log p(y|f,z).
%    This subfunction is needed when using Laplace approximation 
%    or MCMC for inference with non-Gaussian likelihoods. This 
%    subfunction is also used in information criteria (DIC, WAIC) 
%    computations.
%
%  See also
%    LIK_WEIBULL_LLG, LIK_WEIBULL_LLG3, LIK_WEIBULL_LLG2, GPLA_E
  
  if isempty(z)
    error(['lik_weibull -> lik_weibull_ll: missing z!    '... 
           'Weibull likelihood needs the censoring    '...
           'indicators as an extra input z. See, for         '...
           'example, lik_weibull and gpla_e.               ']);
  end

  a = lik.shape;
  ll = sum((1-z).*(log(a) + (a-1).*log(y)-f) - exp(-f).*y.^a);

end

function llg = lik_weibull_llg(lik, y, f, param, z)
%LIK_WEIBULL_LLG  Gradient of the log likelihood
%
%  Description 
%    LLG = LIK_WEIBULL_LLG(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, survival times Y, censoring indicators Z and
%    latent values F. Returns the gradient of the log likelihood
%    with respect to PARAM. At the moment PARAM can be 'param' or
%    'latent'. This subfunction is needed when using Laplace 
%    approximation or MCMC for inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_WEIBULL_LL, LIK_WEIBULL_LLG2, LIK_WEIBULL_LLG3, GPLA_E

  if isempty(z)
    error(['lik_weibull -> lik_weibull_llg: missing z!    '... 
           'Weibull likelihood needs the censoring    '...
           'indicators as an extra input z. See, for         '...
           'example, lik_weibull and gpla_e.               ']);
  end

  a = lik.shape;
  switch param
    case 'param'      
      llg = sum((1-z).*(1./a + log(y)) - exp(-f).*y.^a.*log(y));
      % correction for the log transformation
      llg = llg.*lik.shape;
    case 'latent'
      llg = -(1-z) + exp(-f).*y.^a;
  end
end

function llg2 = lik_weibull_llg2(lik, y, f, param, z)
%LIK_WEIBULL_LLG2  Second gradients of the log likelihood
%
%  Description        
%    LLG2 = LIK_WEIBULL_LLG2(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, survival times Y, censoring indicators Z, and
%    latent values F. Returns the hessian of the log likelihood
%    with respect to PARAM. At the moment PARAM can be only
%    'latent'. LLG2 is a vector with diagonal elements of the
%    Hessian matrix (off diagonals are zero). This subfunction 
%    is needed when using Laplace approximation or EP for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_WEIBULL_LL, LIK_WEIBULL_LLG, LIK_WEIBULL_LLG3, GPLA_E

  if isempty(z)
    error(['lik_weibull -> lik_weibull_llg2: missing z!   '... 
           'Weibull likelihood needs the censoring   '...
           'indicators as an extra input z. See, for         '...
           'example, lik_weibull and gpla_e.               ']);
  end

  a = lik.shape;
  switch param
    case 'param'
      
    case 'latent'
      llg2 = -exp(-f).*y.^a;
    case 'latent+param'
      llg2 = exp(-f).*y.^a.*log(y);
      % correction due to the log transformation
      llg2 = llg2.*lik.shape;
  end
end    

function llg3 = lik_weibull_llg3(lik, y, f, param, z)
%LIK_WEIBULL_LLG3  Third gradients of the log likelihood
%
%  Description
%    LLG3 = LIK_WEIBULL_LLG3(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, survival times Y, censoring indicators Z and
%    latent values F and returns the third gradients of the log
%    likelihood with respect to PARAM. At the moment PARAM can be
%    only 'latent'. LLG3 is a vector with third gradients. This 
%    subfunction is needed when using Laplace approximation for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_WEIBULL_LL, LIK_WEIBULL_LLG, LIK_WEIBULL_LLG2, GPLA_E, GPLA_G

  if isempty(z)
    error(['lik_weibull -> lik_weibull_llg3: missing z!   '... 
           'Weibull likelihood needs the censoring    '...
           'indicators as an extra input z. See, for         '...
           'example, lik_weibull and gpla_e.               ']);
  end

  a = lik.shape;
  switch param
    case 'param'
      
    case 'latent'
      llg3 = exp(-f).*y.^a;
    case 'latent2+param'
      llg3 = -exp(-f).*y.^a.*log(y);
      % correction due to the log transformation
      llg3 = llg3.*lik.shape;
  end
end

function [logM_0, m_1, sigm2hati1] = lik_weibull_tiltedMoments(lik, y, i1, sigm2_i, myy_i, z)
%LIK_WEIBULL_TILTEDMOMENTS  Returns the marginal moments for EP algorithm
%
%  Description
%    [M_0, M_1, M2] = LIK_WEIBULL_TILTEDMOMENTS(LIK, Y, I, S2,
%    MYY, Z) takes a likelihood structure LIK, survival times
%    Y, censoring indicators Z, index I and cavity variance S2 and
%    mean MYY. Returns the zeroth moment M_0, mean M_1 and
%    variance M_2 of the posterior marginal (see Rasmussen and
%    Williams (2006): Gaussian processes for Machine Learning,
%    page 55). This subfunction is needed when using EP for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    GPEP_E
  
%  if isempty(z)
%    error(['lik_weibull -> lik_weibull_tiltedMoments: missing z!'... 
%           'Weibull likelihood needs the censoring            '...
%           'indicators as an extra input z. See, for                 '...
%           'example, lik_weibull and gpep_e.                       ']);
%  end
  
  yy = y(i1);
  yc = 1-z(i1);
  r = lik.shape;
  logM_0=zeros(size(yy));
  m_1=zeros(size(yy));
  sigm2hati1=zeros(size(yy));
  
  for i=1:length(i1)
    % get a function handle of an unnormalized tilted distribution
    % (likelihood * cavity = Weibull * Gaussian)
    % and useful integration limits
    [tf,minf,maxf]=init_weibull_norm(yy(i),myy_i(i),sigm2_i(i),yc(i),r);
    
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
        error('lik_weibull_tilted_moments: sigm2hati1 >= sigm2_i');
      end
    end
    logM_0(i) = log(m_0);
  end
  
end

function [g_i] = lik_weibull_siteDeriv(lik, y, i1, sigm2_i, myy_i, z)
%LIK_WEIBULL_SITEDERIV  Evaluate the expectation of the gradient
%                      of the log likelihood term with respect
%                      to the likelihood parameters for EP 
%
%  Description [M_0, M_1, M2] =
%    LIK_WEIBULL_SITEDERIV(LIK, Y, I, S2, MYY, Z) takes a
%    likelihood structure LIK, survival times Y, expected
%    counts Z, index I and cavity variance S2 and mean MYY. 
%    Returns E_f [d log p(y_i|f_i) /d a], where a is the
%    likelihood parameter and the expectation is over the
%    marginal posterior. This term is needed when evaluating the
%    gradients of the marginal likelihood estimate Z_EP with
%    respect to the likelihood parameters (see Seeger (2008):
%    Expectation propagation for exponential families). This 
%    subfunction is needed when using EP for inference with 
%    non-Gaussian likelihoods and there are likelihood parameters.
%
%  See also
%    GPEP_G

  if isempty(z)
    error(['lik_weibull -> lik_weibull_siteDeriv: missing z!'... 
           'Weibull likelihood needs the censoring        '...
           'indicators as an extra input z. See, for             '...
           'example, lik_weibull and gpla_e.                   ']);
  end

  yy = y(i1);
  yc = 1-z(i1);
  r = lik.shape;
  
  % get a function handle of an unnormalized tilted distribution 
  % (likelihood * cavity = Weibull * Gaussian)
  % and useful integration limits
  [tf,minf,maxf]=init_weibull_norm(yy,myy_i,sigm2_i,yc,r);
  % additionally get function handle for the derivative
  td = @deriv;
  
  % Integrate with quadgk
  [m_0, fhncnt] = quadgk(tf, minf, maxf);
  [g_i, fhncnt] = quadgk(@(f) td(f).*tf(f)./m_0, minf, maxf);
  g_i = g_i.*r;

  function g = deriv(f)
    g = yc.*(1./r + log(yy)) - exp(-f).*yy.^r.*log(yy);
  end
end

function [lpy, Ey, Vary] = lik_weibull_predy(lik, Ef, Varf, yt, zt)
%LIK_WEIBULL_PREDY  Returns the predictive mean, variance and density of y
%
%  Description   
%    LPY = LIK_WEIBULL_PREDY(LIK, EF, VARF YT, ZT)
%    Returns logarithm of the predictive density PY of YT, that is 
%        p(yt | zt) = \int p(yt | f, zt) p(f|y) df.
%    This requires also the survival times YT, censoring indicators ZT.
%    This subfunction is needed when computing posterior predictive 
%    distributions for future observations.
%
%    [LPY, EY, VARY] = LIK_WEIBULL_PREDY(LIK, EF, VARF) takes a
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

  if isempty(zt)
    error(['lik_weibull -> lik_weibull_predy: missing zt!'... 
           'Weibull likelihood needs the censoring    '...
           'indicators as an extra input zt. See, for         '...
           'example, lik_weibull and gpla_e.               ']);
  end

  yc = 1-zt;
  r = lik.shape;
  
  Ey=[];
  Vary=[];
  
  %     lpy = zeros(size(Ef));
  %     Ey = zeros(size(Ef));
  %     EVary = zeros(size(Ef));
  %     VarEy = zeros(size(Ef)); 
  %     
  %     % Evaluate Ey and Vary 
  %     for i1=1:length(Ef)
  %       %%% With quadrature
  %       myy_i = Ef(i1);
  %       sigm_i = sqrt(Varf(i1));
  %       minf=myy_i-6*sigm_i;
  %       maxf=myy_i+6*sigm_i;
  % 
  %       F = @(f) exp(log(yc(i1))+f+norm_lpdf(f,myy_i,sigm_i));
  %       Ey(i1) = quadgk(F,minf,maxf);
  %       
  %       F2 = @(f) exp(log(yc(i1).*exp(f)+((yc(i1).*exp(f)).^2/r))+norm_lpdf(f,myy_i,sigm_i));
  %       EVary(i1) = quadgk(F2,minf,maxf);
  %       
  %       F3 = @(f) exp(2*log(yc(i1))+2*f+norm_lpdf(f,myy_i,sigm_i));
  %       VarEy(i1) = quadgk(F3,minf,maxf) - Ey(i1).^2;
  %     end
  %     Vary = EVary + VarEy;

  % Evaluate the posterior predictive densities of the given observations
  lpy = zeros(length(yt),1);
  for i1=1:length(yt)
    % get a function handle of the likelihood times posterior
    % (likelihood * posterior = Weibull * Gaussian)
    % and useful integration limits
    [pdf,minf,maxf]=init_weibull_norm(...
      yt(i1),Ef(i1),Varf(i1),yc(i1),r);
    % integrate over the f to get posterior predictive distribution
    lpy(i1) = log(quadgk(pdf, minf, maxf));
  end
end

function [df,minf,maxf] = init_weibull_norm(yy,myy_i,sigm2_i,yc,r)
%INIT_WEIBULL_NORM
%
%  Description
%    Return function handle to a function evaluating
%    Weibull * Gaussian which is used for evaluating
%    (likelihood * cavity) or (likelihood * posterior) Return
%    also useful limits for integration. This is private function
%    for lik_weibull. This subfunction is needed by subfunctions
%    tiltedMoments, siteDeriv and predy.
%  
%  See also
%    LIK_WEIBULL_TILTEDMOMENTS, LIK_WEIBULL_SITEDERIV,
%    LIK_WEIBULL_PREDY
  
% avoid repetitive evaluation of constant part
  ldconst = yc*log(r)+yc*(r-1)*log(yy)...
            - log(sigm2_i)/2 - log(2*pi)/2;
  
   
  
  % Create function handle for the function to be integrated
  df = @weibull_norm;
  % use log to avoid underflow, and derivates for faster search
  ld = @log_weibull_norm;
  ldg = @log_weibull_norm_g;
  ldg2 = @log_weibull_norm_g2;

  % Set the limits for integration
  if yc==0
    % with yy==0, the mode of the likelihood is not defined
    % use the mode of the Gaussian (cavity or posterior) as a first guess
    modef = myy_i;
  else
    % use precision weighted mean of the Gaussian approximation
    % of the Weibull likelihood and Gaussian
    mu=-log(yc./(yy.^r));
    %s2=1./(yc+1./sigm2_i);
    s2=1./yc;
    modef = (myy_i/sigm2_i + mu/s2)/(1/sigm2_i + 1/s2);
  end
  % find the mode of the integrand using Newton iterations
  % few iterations is enough, since first guess is in the right direction
  niter=4;       % number of Newton iterations
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
      error(['lik_weibull -> init_weibull_norm: ' ...
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
      error(['lik_weibull -> init_weibull_norm: ' ...
             'integration interval maximun not found ' ...
             'even after looking hard!'])
    end
  end
  
  function integrand = weibull_norm(f)
  % Weibull * Gaussian
    integrand = exp(ldconst ...
                    -yc.*f -exp(-f).*yy.^r ...
                    -0.5*(f-myy_i).^2./sigm2_i);
  end

  function log_int = log_weibull_norm(f)
  % log(Weibull * Gaussian)
  % log_weibull_norm is used to avoid underflow when searching
  % integration interval
    log_int = ldconst ...
              -yc.*f -exp(-f).*yy.^r ...
              -0.5*(f-myy_i).^2./sigm2_i;
  end

  function g = log_weibull_norm_g(f)
  % d/df log(Weibull * Gaussian)
  % derivative of log_weibull_norm
    g = -yc + exp(-f).*yy.^r ...
        + (myy_i - f)./sigm2_i;
  end

  function g2 = log_weibull_norm_g2(f)
  % d^2/df^2 log(Weibull * Gaussian)
  % second derivate of log_weibull_norm
    g2 = - exp(-f).*yy.^r ...
         -1/sigm2_i;
  end

end

function cdf = lik_weibull_predcdf(lik, Ef, Varf, yt)
%LIK_WEIBULL_PREDCDF  Returns the predictive cdf evaluated at yt 
%
%  Description   
%    CDF = LIK_WEIBULL_PREDCDF(LIK, EF, VARF, YT)
%    Returns the predictive cdf evaluated at YT given likelihood
%    structure LIK, posterior mean EF and posterior Variance VARF
%    of the latent variable. This subfunction is needed when using
%    functions gp_predcdf or gp_kfcv_cdf.
%
%  See also
%    GP_PREDCDF

  r = lik.shape;
  
  % Evaluate the posterior predictive cdf at given yt
  cdf = zeros(length(yt),1);
  for i1=1:length(yt)
    % Get a function handle of the likelihood times posterior
    % (likelihood * posterior = Weibull * Gaussian)
    % and useful integration limits.
    % yc=0 when evaluating predictive cdf
    [sf,minf,maxf]=init_weibull_norm(...
      yt(i1),Ef(i1),Varf(i1),0,r);
    % integrate over the f to get posterior predictive distribution
    cdf(i1) = 1-quadgk(sf, minf, maxf);
  end
end

function p = lik_weibull_invlink(lik, f)
%LIK_WEIBULL Returns values of inverse link function
%             
%  Description 
%    P = LIK_WEIBULL_INVLINK(LIK, F) takes a likelihood structure LIK and
%    latent values F and returns the values of inverse link function P.
%    This subfunction is needed when using function gp_predprctmu.
%
%     See also
%     LIK_WEIBULL_LL, LIK_WEIBULL_PREDY

p = exp(f);
end

function reclik = lik_weibull_recappend(reclik, ri, lik)
%RECAPPEND  Append the parameters to the record
%
%  Description 
%    RECLIK = GPCF_WEIBULL_RECAPPEND(RECLIK, RI, LIK) takes a
%    likelihood record structure RECLIK, record index RI and
%    likelihood structure LIK with the current MCMC samples of
%    the parameters. Returns RECLIK which contains all the old
%    samples and the current samples from LIK. This subfunction
%    is needed when using MCMC sampling (gp_mc).
% 
%  See also
%    GP_MC

  if nargin == 2
    % Initialize the record
    reclik.type = 'Weibull';

    % Initialize parameter
    reclik.shape = [];

    % Set the function handles
    reclik.fh.pak = @lik_weibull_pak;
    reclik.fh.unpak = @lik_weibull_unpak;
    reclik.fh.lp = @lik_t_lp;
    reclik.fh.lpg = @lik_t_lpg;
    reclik.fh.ll = @lik_weibull_ll;
    reclik.fh.llg = @lik_weibull_llg;    
    reclik.fh.llg2 = @lik_weibull_llg2;
    reclik.fh.llg3 = @lik_weibull_llg3;
    reclik.fh.tiltedMoments = @lik_weibull_tiltedMoments;
    reclik.fh.invlink = @lik_weibull_invlink;
    reclik.fh.predy = @lik_weibull_predy;
    reclik.fh.predcdf = @lik_weibull_predcdf;
    reclik.fh.recappend = @lik_weibull_recappend;
    reclik.p=[];
    reclik.p.shape=[];
    if ~isempty(ri.p.shape)
      reclik.p.shape = ri.p.shape;
    end
  else
    % Append to the record
    reclik.shape(ri,:)=lik.shape;
    if ~isempty(lik.p)
      reclik.p.shape = lik.p.shape.fh.recappend(reclik.p.shape, ri, lik.p.shape);
    end
  end
end

