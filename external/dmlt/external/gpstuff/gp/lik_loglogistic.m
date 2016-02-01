function lik = lik_loglogistic(varargin)
%LIK_LOGLOGISTIC Create a right censored log-logistic likelihood structure 
%
%  Description
%    LIK = LIK_LOGLOGISTIC('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates a likelihood structure for right censored log-logistic
%    survival model in which the named parameters have the
%    specified values. Any unspecified parameters are set to
%    default values.
%  
%    LIK = LIK_LOGLOGISTIC(LIK,'PARAM1',VALUE1,'PARAM2,VALUE2,...)
%    modify a likelihood structure with the named parameters
%    altered with the specified values.
%
%    Parameters for loggaussian likelihood [default]
%      shape       - shape parameter r [1]
%      shape_prior - prior for shape [prior_logunif]
%  
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc.
%
%    The likelihood is defined as follows:
%                  __ n
%      p(y|f, z) = || i=1 [ (r/exp(f_i)*(y_i/exp(f_i))^(r-1)/
%                               (1+(y_i/exp(f_i))^r))^(1-z_i)
%                               *(1+(y_i/exp(f_i))^r)^(-z_i) ]
%                           
%
%    where r is the shape parameter of log-logistic distribution.
%    z is a vector of censoring indicators with z = 0 for uncensored event
%    and z = 1 for right censored event. 
%
%    When using the loggaussian likelihood you need to give the vector z
%    as an extra parameter to each function that requires also y. 
%    For example, you should call gpla_e as follows: gpla_e(w, gp,
%    x, y, 'z', z)
%
%  See also
%    GP_SET, LIK_*, PRIOR_*
%

% Copyright (c) 2012 Ville Tolvanen

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_LOGLOGISTIC';
  ip.addOptional('lik', [], @isstruct);
  ip.addParamValue('shape',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('shape_prior',prior_logunif(), @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  lik=ip.Results.lik;
  
  if isempty(lik)
    init=true;
    lik.type = 'Log-Logistic';
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'Log-Logistic')
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
    lik.fh.pak = @lik_loglogistic_pak;
    lik.fh.unpak = @lik_loglogistic_unpak;
    lik.fh.lp = @lik_loglogistic_lp;
    lik.fh.lpg = @lik_loglogistic_lpg;
    lik.fh.ll = @lik_loglogistic_ll;
    lik.fh.llg = @lik_loglogistic_llg;    
    lik.fh.llg2 = @lik_loglogistic_llg2;
    lik.fh.llg3 = @lik_loglogistic_llg3;
    lik.fh.tiltedMoments = @lik_loglogistic_tiltedMoments;
    lik.fh.siteDeriv = @lik_loglogistic_siteDeriv;
    lik.fh.invlink = @lik_loglogistic_invlink;
    lik.fh.predy = @lik_loglogistic_predy;
    lik.fh.recappend = @lik_loglogistic_recappend;
    lik.fh.predcdf=@lik_loglogistic_predcdf;
  end

end

function [w,s] = lik_loglogistic_pak(lik)
%LIK_LOGLOGISTIC_PAK  Combine likelihood parameters into one vector.
%
%  Description 
%    W = LIK_LOGLOGISTIC_PAK(LIK) takes a likelihood structure LIK and
%    combines the parameters into a single row vector W. This is a 
%    mandatory subfunction used for example in energy and gradient 
%    computations.
%     
%       w = log(lik.shape)
%
%   See also
%   LIK_LOGLOGISTIC_UNPAK, GP_PAK
  
  w=[];s={};
  if ~isempty(lik.p.shape)
    w = log(lik.shape);
    s = [s; 'log(loglogistic.shape)'];
    [wh sh] = lik.p.shape.fh.pak(lik.p.shape);
    w = [w wh];
    s = [s; sh];
  end
end


function [lik, w] = lik_loglogistic_unpak(lik, w)
%LIK_LOGLOGISTIC_UNPAK  Extract likelihood parameters from the vector.
%
%  Description
%    [LIK, W] = LIK_LOGLOGISTIC_UNPAK(W, LIK) takes a likelihood
%    structure LIK and extracts the parameters from the vector W
%    to the LIK structure. This is a mandatory subfunction used 
%    for example in energy and gradient computations.
%     
%   Assignment is inverse of  
%       w = log(lik.shape)
%
%   See also
%   LIK_LOGLOGISTIC_PAK, GP_UNPAK

  if ~isempty(lik.p.shape)
    lik.shape = exp(w(1));
    w = w(2:end);
    [p, w] = lik.p.shape.fh.unpak(lik.p.shape, w);
    lik.p.shape = p;
  end
end


function lp = lik_loglogistic_lp(lik, varargin)
%LIK_LOGLOGISTIC_LP  log(prior) of the likelihood parameters
%
%  Description
%    LP = LIK_LOGLOGISTIC_LP(LIK) takes a likelihood structure LIK and
%    returns log(p(th)), where th collects the parameters. This subfunction
%    is needed when there are likelihood parameters.
%
%  See also
%    LIK_LOGLOGISTIC_LLG, LIK_LOGLOGISTIC_LLG3, LIK_LOGLOGISTIC_LLG2, GPLA_E
  

% If prior for shape parameter, add its contribution
  lp=0;
  if ~isempty(lik.p.shape)
    lp = lik.p.shape.fh.lp(lik.shape, lik.p.shape) +log(lik.shape);
  end
  
end


function lpg = lik_loglogistic_lpg(lik)
%LIK_LOGLOGISTIC_LPG  d log(prior)/dth of the likelihood 
%                parameters th
%
%  Description
%    E = LIK_LOGLOGISTIC_LPG(LIK) takes a likelihood structure LIK and
%    returns d log(p(th))/dth, where th collects the parameters. This
%    subfunction is needed when there are likelihood parameters.
%
%  See also
%    LIK_LOGLOGISTIC_LLG, LIK_LOGLOGISTIC_LLG3, LIK_LOGLOGISTIC_LLG2, GPLA_G
  
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

function ll = lik_loglogistic_ll(lik, y, f, z)
%LIK_LOGLOGISTIC_LL  Log likelihood
%
%  Description
%    LL = LIK_LOGLOGISTIC_LL(LIK, Y, F, Z) takes a likelihood
%    structure LIK, survival times Y, censoring indicators Z, and
%    latent values F. Returns the log likelihood, log p(y|f,z).
%    This subfunction is needed when using Laplace approximation 
%    or MCMC for inference with non-Gaussian likelihoods. This 
%    subfunction is also used in information criteria (DIC, WAIC) 
%    computations.
%
%  See also
%    LIK_LOGLOGISTIC_LLG, LIK_LOGLOGISTIC_LLG3, LIK_LOGLOGISTIC_LLG2, GPLA_E
  
  if isempty(z)
    error(['lik_loglogistic -> lik_loglogistic_ll: missing z!    '... 
           'loglogistic likelihood needs the censoring    '...
           'indicators as an extra input z. See, for         '...
           'example, lik_loglogistic and gpla_e.               ']);
  end

  r = lik.shape;
  ll = sum((1-z).*(log(r)+(r-1).*log(y)-r.*f) +(z-2).*log(1+(y./exp(f)).^r));

end

function llg = lik_loglogistic_llg(lik, y, f, param, z)
%LIK_LOGLOGISTIC_LLG  Gradient of the log likelihood
%
%  Description 
%    LLG = LIK_LOGLOGISTIC_LLG(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, survival times Y, censoring indicators Z and
%    latent values F. Returns the gradient of the log likelihood
%    with respect to PARAM. At the moment PARAM can be 'param' or
%    'latent'. This subfunction is needed when using Laplace 
%    approximation or MCMC for inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_LOGLOGISTIC_LL, LIK_LOGLOGISTIC_LLG2, LIK_LOGLOGISTIC_LLG3, GPLA_E

  if isempty(z)
    error(['lik_loglogistic -> lik_loglogistic_llg: missing z!    '... 
           'loglogistic likelihood needs the censoring    '...
           'indicators as an extra input z. See, for         '...
           'example, lik_loglogistic and gpla_e.               ']);
  end

  r = lik.shape;
  switch param
    case 'param'      
      llg = sum((1-z).*(1/r+log(y)-f) + (z-2)./(1+(y./exp(f)).^r).* ...
             (y./exp(f)).^r.*(log(y)-f));
      % correction for the log transformation
      llg = llg.*lik.shape;
    case 'latent'
      llg = -r.*(1-z) - (z-2).*r.*(y./exp(f)).^r./(1+(y./exp(f)).^r);
  end
end

function llg2 = lik_loglogistic_llg2(lik, y, f, param, z)
%LIK_LOGLOGISTIC_LLG2  Second gradients of the log likelihood
%
%  Description        
%    LLG2 = LIK_LOGLOGISTIC_LLG2(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, survival times Y, censoring indicators Z, and
%    latent values F. Returns the hessian of the log likelihood
%    with respect to PARAM. At the moment PARAM can be only
%    'latent'. LLG2 is a vector with diagonal elements of the
%    Hessian matrix (off diagonals are zero). This subfunction 
%    is needed when using Laplace approximation or EP for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_LOGLOGISTIC_LL, LIK_LOGLOGISTIC_LLG, LIK_LOGLOGISTIC_LLG3, GPLA_E

  if isempty(z)
    error(['lik_loglogistic -> lik_loglogistic_llg2: missing z!   '... 
           'loglogistic likelihood needs the censoring   '...
           'indicators as an extra input z. See, for         '...
           'example, lik_loglogistic and gpla_e.               ']);
  end

  r = lik.shape;
  switch param
    case 'param'
      
    case 'latent'
      llg2 = r.^2.*(z-2).*(y./exp(f)).^r./(1+(y./exp(f)).^r).^2;
    case 'latent+param'
      llg2 = (z-1) - (z-2).*(y./exp(f)).^r./(1+(y./exp(f)).^r) ...
               + (z-2).*r.*(y./exp(f)).^(2*r).*(log(y)-f)./(1+(y./exp(f)).^r).^2 ...
               - (z-2).*r.*(y./exp(f)).^r.*(log(y)-f)./(1+(y./exp(f)).^r);
      % correction due to the log transformation
      llg2 = llg2.*r;
  end
end    

function llg3 = lik_loglogistic_llg3(lik, y, f, param, z)
%LIK_LOGLOGISTIC_LLG3  Third gradients of the log likelihood
%
%  Description
%    LLG3 = LIK_LOGLOGISTIC_LLG3(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, survival times Y, censoring indicators Z and
%    latent values F and returns the third gradients of the log
%    likelihood with respect to PARAM. At the moment PARAM can be
%    only 'latent'. LLG3 is a vector with third gradients. This 
%    subfunction is needed when using Laplace approximation for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_LOGLOGISTIC_LL, LIK_LOGLOGISTIC_LLG, LIK_LOGLOGISTIC_LLG2, GPLA_E, GPLA_G

  if isempty(z)
    error(['lik_loglogistic -> lik_loglogistic_llg3: missing z!   '... 
           'loglogistic likelihood needs the censoring    '...
           'indicators as an extra input z. See, for         '...
           'example, lik_loglogistic and gpla_e.               ']);
  end

  r = lik.shape;
  switch param
    case 'param'
      
    case 'latent'
      llg3 = r.^3.*(z-2).*(y./exp(f)).^r.*(-1+(y./exp(f)).^r)./(1+(y./exp(f)).^r).^3;
    case 'latent2+param'
      llg3 = -(r.*(z-2).*(y./exp(f)).^r.*(-2-2.*(y./exp(f)).^r + ...
              r.*(-1+(y./exp(f)).^r).*log(y./exp(f))))./(1+(y./exp(f)).^r).^3;
      % correction due to the log transformation
      llg3 = llg3.*lik.shape;
  end
end

function [logM_0, m_1, sigm2hati1] = lik_loglogistic_tiltedMoments(lik, y, i1, sigm2_i, myy_i, z)
%LIK_LOGLOGISTIC_TILTEDMOMENTS  Returns the marginal moments for EP algorithm
%
%  Description
%    [M_0, M_1, M2] = LIK_LOGLOGISTIC_TILTEDMOMENTS(LIK, Y, I, S2,
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
  
 if isempty(z)
   error(['lik_loglogistic -> lik_loglogistic_tiltedMoments: missing z!'... 
          'loglogistic likelihood needs the censoring            '...
          'indicators as an extra input z. See, for                 '...
          'example, lik_loglogistic and gpep_e.                       ']);
 end
  
  yy = y(i1);
  yc = 1-z(i1);
  r = lik.shape;
  logM_0=zeros(size(yy));
  m_1=zeros(size(yy));
  sigm2hati1=zeros(size(yy));
  
  for i=1:length(i1)
    % get a function handle of an unnormalized tilted distribution
    % (likelihood * cavity = Negative-binomial * Gaussian)
    % and useful integration limits
    [tf,minf,maxf]=init_loglogistic_norm(yy(i),myy_i(i),sigm2_i(i),yc(i),r);
    
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
        error('lik_loglogistic_tilted_moments: sigm2hati1 >= sigm2_i');
      end
    end
    logM_0(i) = log(m_0);
  end
end

function [g_i] = lik_loglogistic_siteDeriv(lik, y, i1, sigm2_i, myy_i, z)
%LIK_LOGLOGISTIC_SITEDERIV  Evaluate the expectation of the gradient
%                      of the log likelihood term with respect
%                      to the likelihood parameters for EP 
%
%  Description [M_0, M_1, M2] =
%    LIK_LOGLOGISTIC_SITEDERIV(LIK, Y, I, S2, MYY, Z) takes a
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
    error(['lik_loglogistic -> lik_loglogistic_siteDeriv: missing z!'... 
           'loglogistic likelihood needs the censoring        '...
           'indicators as an extra input z. See, for             '...
           'example, lik_loglogistic and gpla_e.                   ']);
  end

  yy = y(i1);
  yc = 1-z(i1);
  r = lik.shape;
  
  % get a function handle of an unnormalized tilted distribution 
  % (likelihood * cavity = Log-Gaussian * Gaussian)
  % and useful integration limits
  [tf,minf,maxf]=init_loglogistic_norm(yy,myy_i,sigm2_i,yc,r);
  % additionally get function handle for the derivative
  td = @deriv;
  
  % Integrate with quadgk
  [m_0, fhncnt] = quadgk(tf, minf, maxf);
  [g_i, fhncnt] = quadgk(@(f) td(f).*tf(f)./m_0, minf, maxf);
  g_i = g_i.*r;

  function g = deriv(f)
    g = yc.*(1/r+log(yy)-f) + (-1-yc)./(1+(yy./exp(f)).^r).* ...
             (yy./exp(f)).^r.*(log(yy)-f);
  end
end

function [lpy, Ey, Vary] = lik_loglogistic_predy(lik, Ef, Varf, yt, zt)
%LIK_LOGLOGISTIC_PREDY  Returns the predictive mean, variance and density of y
%
%  Description   
%    LPY = LIK_LOGLOGISTIC_PREDY(LIK, EF, VARF YT, ZT)
%    Returns logarithm of the predictive density PY of YT, that is 
%        p(yt | zt) = \int p(yt | f, zt) p(f|y) df.
%    This requires also the survival times YT, censoring indicators ZT.
%    This subfunction is needed when computing posterior predictive 
%    distributions for future observations.
%
%    [LPY, EY, VARY] = LIK_LOGLOGISTIC_PREDY(LIK, EF, VARF) takes a
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
    error(['lik_loglogistic -> lik_loglogistic_predy: missing zt!'... 
           'loglogistic likelihood needs the censoring    '...
           'indicators as an extra input zt. See, for         '...
           'example, lik_loglogistic and gpla_e.               ']);
  end

  yc = 1-zt;
  r = lik.shape;
  
  Ey=[];
  Vary=[];

  % Evaluate the posterior predictive densities of the given observations
  lpy = zeros(length(yt),1);
  for i1=1:length(yt)
    if abs(Ef(i1))>700
      lpy(i1) = NaN;
    else
      % get a function handle of the likelihood times posterior
      % (likelihood * posterior = Negative-binomial * Gaussian)
      % and useful integration limits
      [pdf,minf,maxf]=init_loglogistic_norm(...
        yt(i1),Ef(i1),Varf(i1),yc(i1),r);
      % integrate over the f to get posterior predictive distribution
      lpy(i1) = log(quadgk(pdf, minf, maxf));
    end
  end
end

function [df,minf,maxf] = init_loglogistic_norm(yy,myy_i,sigm2_i,yc,r)
%INIT_LOGLOGISTIC_NORM
%
%  Description
%    Return function handle to a function evaluating
%    loglogistic * Gaussian which is used for evaluating
%    (likelihood * cavity) or (likelihood * posterior) Return
%    also useful limits for integration. This is private function
%    for lik_loglogistic. This subfunction is needed by subfunctions
%    tiltedMoments, siteDeriv and predy.
%  
%  See also
%    LIK_LOGLOGISTIC_TILTEDMOMENTS, LIK_LOGLOGISTIC_SITEDERIV,
%    LIK_LOGLOGISTIC_PREDY
  
% avoid repetitive evaluation of constant part
  ldconst =  yc.*log(r) + yc.*(r-1).*log(yy) ...
            - log(sigm2_i)/2 - log(2*pi)/2;
  
  % Create function handle for the function to be integrated
  df = @loglogistic_norm;
  % use log to avoid underflow, and derivates for faster search
  ld = @log_loglogistic_norm;
  ldg = @log_loglogistic_norm_g;
  ldg2 = @log_loglogistic_norm_g2;

  % Set the limits for integration
  if yc==0
    % with yy==0, the mode of the likelihood is not defined
    % use the mode of the Gaussian (cavity or posterior) as a first guess
    modef = myy_i;
  else
    % use precision weighted mean of the Gaussian approximation
    % of the loglogistic likelihood and Gaussian
    mu=-log(yc.^(1/r)./yy);
    s2=-r^2.*yc.*(-1-yc)./(1+yc).^2;
%     s2=1;
    modef = (myy_i/sigm2_i + mu/s2)/(1/sigm2_i + 1/s2);
  end
  % find the mode of the integrand using Newton iterations
  % few iterations is enough, since the first guess in the right direction
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
      error(['lik_loglogistic -> init_loglogistic_norm: ' ...
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
      error(['lik_loglogistic -> init_loglogistic_norm: ' ...
             'integration interval maximun not found ' ...
             'even after looking hard!'])
    end
  end
  
  function integrand = loglogistic_norm(f)
  % loglogistic * Gaussian
    integrand = exp(ldconst ...
                     - yc.*r.*f +(-1-yc).*log(1+(yy./exp(f)).^r) ...
                    -0.5*(f-myy_i).^2./sigm2_i);
  end

  function log_int = log_loglogistic_norm(f)
  % log(loglogistic * Gaussian)
  % log_loglogistic_norm is used to avoid underflow when searching
  % integration interval
    log_int = ldconst ...
               -yc.*r.*f +(-1-yc).*log(1+(yy./exp(f)).^r) ...
              -0.5*(f-myy_i).^2./sigm2_i;
  end

  function g = log_loglogistic_norm_g(f)
  % d/df log(loglogistic * Gaussian)
  % derivative of log_loglogistic_norm
    g = -r.*yc - (-1-yc).*r.*(yy./exp(f)).^r./(1+(yy./exp(f)).^r) ...
        + (myy_i - f)./sigm2_i;
  end

  function g2 = log_loglogistic_norm_g2(f)
  % d^2/df^2 log(loglogistic * Gaussian)
  % second derivate of log_loglogistic_norm
    g2 =  r.^2.*(-1-yc).*(yy./exp(f)).^r./(1+(yy./exp(f)).^r).^2 ...
              -1/sigm2_i;
  end

end

function cdf = lik_loglogistic_predcdf(lik, Ef, Varf, yt)
%LIK_LOGLOGISTIC_PREDCDF  Returns the predictive cdf evaluated at yt 
%
%  Description   
%    CDF = LIK_LOGLOGISTIC_PREDCDF(LIK, EF, VARF, YT)
%    Returns the predictive cdf evaluated at YT given likelihood
%    structure LIK, posterior mean EF and posterior Variance VARF
%    of the latent variable. This subfunction is needed when using
%    functions gp_predcdf or gp_kfcv_cdf.
%
%  See also
%    GP_PREDCDF

  r = lik.shape;
  
  % Evaluate the posterior predictive densities of the given observations
  cdf = zeros(length(yt),1);
  for i1=1:length(yt)
    % Get a function handle of the likelihood times posterior
    % (likelihood * posterior = log-logistic * Gaussian)
    % and useful integration limits.
    % yc=0 when evaluating predictive cdf
    [pdf,minf,maxf]=init_loglogistic_norm(...
      yt(i1),Ef(i1),Varf(i1),0,r);
    % integrate over the f to get posterior predictive distribution
    cdf(i1) = 1-quadgk(pdf, minf, maxf);
  end
end

function p = lik_loglogistic_invlink(lik, f)
%LIK_loglogistic Returns values of inverse link function
%             
%  Description 
%    P = LIK_LOGLOGISTIC_INVLINK(LIK, F) takes a likelihood structure LIK and
%    latent values F and returns the values of inverse link function P.
%    This subfunction is needed when using function gp_predprctmu.
%
%     See also
%     LIK_LOGLOGISTIC_LL, LIK_LOGLOGISTIC_PREDY

p = exp(f);
end

function reclik = lik_loglogistic_recappend(reclik, ri, lik)
%RECAPPEND  Append the parameters to the record
%
%  Description 
%    RECLIK = GPCF_LOGLOGISTIC_RECAPPEND(RECLIK, RI, LIK) takes a
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
    reclik.type = 'Log-Logistic';

    % Initialize parameter
    reclik.shape = [];

    % Set the function handles
    reclik.fh.pak = @lik_loglogistic_pak;
    reclik.fh.unpak = @lik_loglogistic_unpak;
    reclik.fh.lp = @lik_loglogistic_lp;
    reclik.fh.lpg = @lik_loglogistic_lpg;
    reclik.fh.ll = @lik_loglogistic_ll;
    reclik.fh.llg = @lik_loglogistic_llg;    
    reclik.fh.llg2 = @lik_loglogistic_llg2;
    reclik.fh.llg3 = @lik_loglogistic_llg3;
    reclik.fh.tiltedMoments = @lik_loglogistic_tiltedMoments;
    reclik.fh.invlink = @lik_loglogistic_invlink;
    reclik.fh.predy = @lik_loglogistic_predy;
    reclik.fh.predcdf=@lik_loglogistic_predcdf;
    reclik.fh.recappend = @lik_loglogistic_recappend;
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
