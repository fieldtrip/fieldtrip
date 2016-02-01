function lik = lik_binomial(varargin)
%LIK_BINOMIAL  Create a Binomial likelihood structure 
%
%  Description
%    LIK = LIK_BINOMIAL creates Binomial likelihood structure.
%
%    The likelihood is defined as follows:
%                  __ n
%      p(y|f, z) = || i=1 [ p_i^(y_i)*(1-p_i)^(z_i-y_i)) * 
%                           gamma(z_i+1)/(gamma(y_i+1)*gamma(z_i-y_i+1))]
%    where p_i = exp(f_i)/ (1+exp(f_i)) is the succes probability,
%    which is a function of the latent variable f_i and z is a
%    vector of numbers of trials. 
%
%    When using Binomial likelihood you need to give the vector z
%    as an extra parameter to each function that requires y also. 
%    For example, you should call gpla_e as follows
%      gpla_e(w, gp, x, y, 'z', z)
%
%  See also
%    GP_SET, LIK_*
%

% Copyright (c) 2009-2010 Jaakko Riihim�ki & Jarno Vanhatalo
% Copyright (c) 2010-2011 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_BINOMIAL';
  ip.addOptional('lik', [], @isstruct);
  ip.parse(varargin{:});
  lik=ip.Results.lik;

  if isempty(lik)
    init=true;
    lik.type = 'Binomial';
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'Binomial')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end

  if init
    % Set the function handles to the subfunctions
    lik.fh.pak = @lik_binomial_pak;
    lik.fh.unpak = @lik_binomial_unpak;
    lik.fh.ll = @lik_binomial_ll;
    lik.fh.llg = @lik_binomial_llg;    
    lik.fh.llg2 = @lik_binomial_llg2;
    lik.fh.llg3 = @lik_binomial_llg3;
    lik.fh.tiltedMoments = @lik_binomial_tiltedMoments;
    lik.fh.predy = @lik_binomial_predy;
    lik.fh.predprcty = @lik_binomial_predprcty;
    lik.fh.invlink = @lik_binomial_invlink;
    lik.fh.recappend = @lik_binomial_recappend;
  end

end

function [w,s] = lik_binomial_pak(lik)
%LIK_BINOMIAL_PAK  Combine likelihood parameters into one vector.
%
%  Description 
%    W = LIK_BINOMIAL_PAK(LIK) takes a likelihood structure LIK
%    and returns an empty verctor W. If Binomial likelihood had
%    parameters this would combine them into a single row vector
%    W (see e.g. likelih_negbin). This is a mandatory subfunction 
%    used for example in energy and gradient computations.
%
%  See also
%    LIK_NEGBIN_UNPAK, GP_PAK

  w = []; s = {};
end


function [lik, w] = lik_binomial_unpak(lik, w)
%LIK_BINOMIAL_UNPAK  Extract likelihood parameters from the vector.
%
%  Description
%    W = LIK_BINOMIAL_UNPAK(W, LIK) Doesn't do anything.
% 
%    If Binomial likelihood had parameters this would extracts
%    them parameters from the vector W to the LIK structure. 
%    This is a mandatory subfunction used for example in energy 
%    and gradient computations.
%
%  See also
%    LIK_BINOMIAL_PAK, GP_UNPAK

  lik=lik;
  w=w;
  
end



function ll = lik_binomial_ll(lik, y, f, z)
%LIK_BINOMIAL_LL  Log likelihood
%
%  Description
%    LL = LIK_BINOMIAL_LL(LIK, Y, F, Z) takes a likelihood
%    structure LIK, succes counts Y, numbers of trials Z, and
%    latent values F. Returns the log likelihood, log p(y|f,z).
%    This subfunction is needed when using Laplace approximation
%    or MCMC for inference with non-Gaussian likelihoods. This 
%    subfunction is also used in information criteria (DIC, WAIC)
%    computations.
%
%  See also
%    LIK_BINOMIAL_LLG, LIK_BINOMIAL_LLG3, LIK_BINOMIAL_LLG2, GPLA_E
  
  if isempty(z)
    error(['lik_binomial -> lik_binomial_ll: missing z!'... 
           'Binomial likelihood needs the expected number of   '...
           'occurrences as an extra input z. See, for         '...
           'example, lik_binomial and gpla_e.             ']);
  end
  
  expf = exp(f);
  p = expf ./ (1+expf);
  N = z;
  ll =  sum(gammaln(N+1)-gammaln(y+1)-gammaln(N-y+1)+y.*log(p)+(N-y).*log(1-p));
end


function llg = lik_binomial_llg(lik, y, f, param, z)
%LIK_BINOMIAL_LLG    Gradient of the log likelihood
%
%  Description 
%    LLG = LIK_BINOMIAL_LLG(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, succes counts Y, numbers of trials Z and
%    latent values F. Returns the gradient of the log likelihood
%    with respect to PARAM. At the moment PARAM can be 'param' or
%    'latent'. This subfunction is needed when using Laplace 
%    approximation or MCMC for inference with non-Gaussian 
%    likelihoods.
%
%  See also
%    LIK_BINOMIAL_LL, LIK_BINOMIAL_LLG2, LIK_BINOMIAL_LLG3, GPLA_E

  if isempty(z)
    error(['lik_binomial -> lik_binomial_llg: missing z!'... 
           'Binomial likelihood needs the expected number of   '...
           'occurrences as an extra input z. See, for         '...
           'example, lik_binomial and gpla_e.             ']);
  end
  
  switch param
    case 'latent'
      expf = exp(f);
      N = z;
      
      llg = y./(1+expf) - (N-y).*expf./(1+expf);
  end
end


function llg2 = lik_binomial_llg2(lik, y, f, param, z)
%LIK_BINOMIAL_LLG2  Second gradients of the log likelihood
%
%  Description        
%    LLG2 = LIK_BINOMIAL_LLG2(LIK, Y, F, PARAM) takes a
%    likelihood structure LIK, succes counts Y, numbers of trials
%    Z, and latent values F. Returns the Hessian of the log
%    likelihood with respect to PARAM. At the moment PARAM can be
%    only 'latent'. G2 is a vector with diagonal elements of the
%    Hessian matrix (off diagonals are zero). This subfunction
%    is needed when using Laplace approximation or EP for inference 
%    with non-Gaussian likelihoods.
%
%  See also
%    LIK_BINOMIAL_LL, LIK_BINOMIAL_LLG, LIK_BINOMIAL_LLG3, GPLA_E

  if isempty(z)
    error(['lik_binomial -> lik_binomial_llg2: missing z!'... 
           'Binomial likelihood needs the expected number of    '...
           'occurrences as an extra input z. See, for          '...
           'example, lik_binomial and gpla_e.              ']);
  end
  
  switch param
    case 'latent'
      expf = exp(f);
      N = z;

      llg2 = -N.*expf./(1+expf).^2;
  end
end


function llg3 = lik_binomial_llg3(lik, y, f, param, z)
%LIK_BINOMIAL_LLG3  Third gradients of the log likelihood
%
%  Description
%    LLG3 = LIK_BINOMIAL_LLG3(LIK, Y, F, PARAM) takes a
%    likelihood structure LIK, succes counts Y, numbers of trials
%    Z and latent values F and returns the third gradients of the
%    log likelihood with respect to PARAM. At the moment PARAM
%    can be only 'latent'. G3 is a vector with third gradients.
%    This subfunction is needed when using Laplace appoximation 
%    for inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_BINOMIAL_LL, LIK_BINOMIAL_LLG, LIK_BINOMIAL_LLG2, GPLA_E, GPLA_G
  
  if isempty(z)
    error(['lik_binomial -> lik_binomial_llg3: missing z!'... 
           'Binomial likelihood needs the expected number of    '...
           'occurrences as an extra input z. See, for          '...
           'example, lik_binomial and gpla_e.              ']);
  end
  
  switch param
    case 'latent'
      expf = exp(f);
      N = z;
      llg3 = N.*(expf.*(expf-1))./(1+expf).^3;
  end
end

function [logM_0, m_1, sigm2hati1] = lik_binomial_tiltedMoments(lik, y, i1, sigm2_i, myy_i, z)
%LIK_BINOMIAL_TILTEDMOMENTS  Returns the marginal moments for EP algorithm
%
%  Description
%    [M_0, M_1, M2] = LIK_BINOMIAL_TILTEDMOMENTS(LIK, Y, I, S2,
%    MYY, Z) takes a likelihood structure LIK, succes counts Y,
%    numbers of trials Z, index I and cavity variance S2 and mean
%    MYY. Returns the zeroth moment M_0, mean M_1 and variance
%    M_2 of the posterior marginal (see Rasmussen and Williams
%    (2006): Gaussian processes for Machine Learning, page 55).
%    This subfunction is needed when using EP for inference with
%    non-Gaussian likelihoods.
%
%  See also
%    GPEP_E
  
%  if isempty(z)
%    error(['lik_binomial -> lik_binomial_tiltedMoments: missing z!'... 
%           'Binomial likelihood needs the expected number of               '...
%           'occurrences as an extra input z. See, for                     '...
%           'example, lik_binomial and gpla_e.                         ']);
%  end
  
  yy = y(i1);
  N = z(i1);
  logM_0=zeros(size(yy));
  m_1=zeros(size(yy));
  sigm2hati1=zeros(size(yy));  
  
  for i=1:length(i1)
    % Create function handle for the function to be integrated
    % (likelihood * cavity) and useful integration limits
    [tf,minf,maxf]=init_binomial_norm(yy(i),myy_i(i),sigm2_i(i),N(i));
    
    % Integrate with quadrature
    RTOL = 1.e-6;
    ATOL = 1.e-10;
    [m_0, m_1(i), m_2] = quad_moments(tf,minf, maxf, RTOL, ATOL);
    sigm2hati1(i) = m_2 - m_1(i).^2;
    
    % If the second central moment is less than cavity variance
    % integrate more precisely. Theoretically for log-concave
    % likelihood should be sigm2hati1 < sigm2_i.
    if sigm2hati1(i) >= sigm2_i(i)
      ATOL = ATOL.^2;
      RTOL = RTOL.^2;
      [m_0, m_1(i), m_2] = quad_moments(tf, minf, maxf, RTOL, ATOL);
      sigm2hati1 = m_2 - m_1(i).^2;
      %    if sigm2hati1 >= sigm2_i
      %      error('lik_binomial_tilted_moments: sigm2hati1 >= sigm2_i');
      %    end
    end
    logM_0(i) = log(m_0);
  end
end


function [lpy, Ey, Vary] = lik_binomial_predy(lik, Ef, Varf, yt, zt)
%LIK_BINOMIAL_PREDY  Returns the predictive mean, variance and density of y
%
%  Description         
%    [LPY] = LIK_BINOMIAL_PREDY(LIK, EF, VARF YT, ZT)
%    Returns logarithm of the predictive density PY of YT, that is 
%        p(yt | y, zt) = \int p(yt | f, zt) p(f|y) df.
%    This requires also the succes counts YT, numbers of trials ZT.
%    This subfunction is needed when computing posterior predictive 
%    distributions for future observations.
%
%    [LPY, EY, VARY] = LIK_BINOMIAL_PREDY(LIK, EF, VARF) takes a
%    likelihood structure LIK, posterior mean EF and posterior
%    Variance VARF of the latent variable and returns the
%    posterior predictive mean EY and variance VARY of the
%    observations related to the latent variables. This subfunction 
%    is needed when computing posterior predictive distributions for 
%    future observations.
%        
%
%  See also 
%    GPEP_PRED, GPLA_PRED, GPMC_PRED

  if isempty(zt)
    error(['lik_binomial -> lik_binomial_predy: missing z!'... 
           'Binomial likelihood needs the expected number of       '...
           'occurrences as an extra input z. See, for             '...
           'example, lik_binomial and gpla_e.                 ']);
  end
  
  if nargout > 1
    nt=length(Ef);
    Ey=zeros(nt,1);
    EVary = zeros(nt,1);
    VarEy = zeros(nt,1);
    for i1=1:nt
      ci = sqrt(Varf(i1));
      F  = @(x)zt(i1)./(1+exp(-x)).*norm_pdf(x,Ef(i1),sqrt(Varf(i1)));
      Ey(i1) = quadgk(F,Ef(i1)-6*ci,Ef(i1)+6*ci);
      
      F2  = @(x)zt(i1)./(1+exp(-x)).*(1-1./(1+exp(-x))).*norm_pdf(x,Ef(i1),sqrt(Varf(i1)));
      EVary(i1) = quadgk(F2,Ef(i1)-6*ci,Ef(i1)+6*ci);
      
      F3  = @(x)(zt(i1)./(1+exp(-x))).^2.*norm_pdf(x,Ef(i1),sqrt(Varf(i1)));
      VarEy(i1) = quadgk(F3,Ef(i1)-6*ci,Ef(i1)+6*ci) - Ey(i1).^2;
    end
    Vary = EVary+VarEy;
  end
  
  nt=length(yt);
  lpy=zeros(nt,1);
  for i1=1:nt
    ci = sqrt(Varf(i1));
    F  = @(x)exp(gammaln(zt(i1)+1)-gammaln(yt(i1)+1)-gammaln(zt(i1)-yt(i1)+1) + yt(i1).*log(1./(1+exp(-x))) + (zt(i1)-yt(i1)).*log(1-(1./(1+exp(-x))))).*norm_pdf(x,Ef(i1),sqrt(Varf(i1)));
    lpy(i1) = log(quadgk(F,Ef(i1)-6*ci,Ef(i1)+6*ci));
  end
  
end

function prctys = lik_binomial_predprcty(lik, Ef, Varf, zt, prcty)
%LIK_BINOMIAL_PREDPRCTY  Returns the percentiled of predictive density of y
%
%  Description         
%    PRCTY = LIK_BINOMIAL_PREDPRCTY(LIK, EF, VARF YT, ZT)
%    Returns percentiles of the predictive density PY of YT, that is 
%    This requires also the succes counts YT, numbers of trials ZT. This
%    subfunction is needed when using function gp_predprcty.
%
%  See also 
%    GP_PREDPCTY

  if isempty(zt)
    error(['lik_binomial -> lik_binomial_predprcty: missing z!'... 
           'Binomial likelihood needs the expected number of       '...
           'occurrences as an extra input z. See, for             '...
           'example, lik_binomial and gpla_e.                 ']);
  end
  
  opt=optimset('TolX',.5,'Display','off');
  nt=size(Ef,1);
  prctys = zeros(nt,numel(prcty));
  prcty=prcty/100;
  for i1=1:nt
    ci = sqrt(Varf(i1));
    for i2=1:numel(prcty)
      a=floor(fminbnd(@(a) (quadgk(@(f) binocdf(a,zt(i1),logitinv(f)).*norm_pdf(f,Ef(i1),ci),Ef(i1)-6*ci,Ef(i1)+6*ci,'AbsTol',1e-4)-prcty(i2)).^2,binoinv(prcty(i2),zt(i1),logitinv(Ef(i1)-1.96*ci)),binoinv(prcty(i2),zt(i1),logitinv(Ef(i1)+1.96*ci)),opt));
      if quadgk(@(f) binocdf(a,zt(i1),logitinv(f)).*norm_pdf(f,Ef(i1),ci),Ef(i1)-6*ci,Ef(i1)+6*ci,'AbsTol',1e-4)<prcty(i2)
        a=a+1;
      end
      prctys(i1,i2)=a;
    end
  end
end

function [df,minf,maxf] = init_binomial_norm(yy,myy_i,sigm2_i,N)
%INIT_LOGIT_NORM
%
%  Description
%    Return function handle to a function evaluating Binomial *
%    Gaussian which is used for evaluating (likelihood * cavity)
%    or (likelihood * posterior) Return also useful limits for
%    integration. This is private function for lik_binomial. This
%    subfunction is needed by subfunctions tiltedMoments and predy.
%  
% See also
%   LIK_BINOMIAL_TILTEDMOMENTS, LIK_BINOMIAL_PREDY
  
% avoid repetitive evaluation of constant part
  ldconst = gammaln(N+1)-gammaln(yy+1)-gammaln(N-yy+1) - log(sigm2_i)/2 - log(2*pi)/2;
%   ldconst = log(factorial(N)/(factorial(yy)*factorial(N-yy))-log(sigm2_i)/2 -log(2*pi)/2;
  
 % Create function handle for the function to be integrated
  df = @binomial_norm;
 % use log to avoid underflow, and derivates for faster search
  ld = @log_binomial_norm;
  ldg = @log_binomial_norm_g;
  ldg2 = @log_binomial_norm_g2;
  
  % Set the limits for integration
  % Binomial likelihood is log-concave so the binomial_norm
  % function is unimodal, which makes things easier
  if yy==0 || yy==N
    % with yy==0 or yy==N the mode of the likelihood is not defined
    % use the mode of the Gaussian (cavity or posterior) as a first guess
    modef = myy_i;
  else
    % use precision weighted mean of the Gaussian approximation of the
    % binomial likelihood and Gaussian
    mean_app = log(yy./(N-yy));
    ld0=1/(1+exp(-mean_app));
    ld1=(1-ld0)*ld0;
    ld2=ld0-3*ld0^2+2*ld0^3;
    var_app=inv(-( yy*(ld2*ld0-ld1^2)/ld0^2 + (N-yy)*(ld2*(ld0-1)-ld1^2)/(ld0-1)^2 ));
    
    modef = (myy_i/sigm2_i + mean_app/var_app)/(1/sigm2_i + 1/var_app);
%     sigm_app = sqrt((1/sigm2_i + 1/var_app)^-1);
  end
  % find the mode of the integrand using Newton iterations
  % few iterations is enough, since the first guess in the right direction
  niter=3;       % number of Newton iterations
  mindelta=1e-6; % tolerance in stopping Newton iterations
  for ni=1:niter
      g = ldg(modef);
      h = ldg2(modef);
      delta=-g/h;
      modef=modef+delta;
      if abs(delta)<mindelta
          break
      end
  end
  % integrand limits based on Gaussian approximation at mode
  modes=sqrt(-1/h);
  minf=modef-4*modes;
  maxf=modef+4*modes;
  modeld=ld(modef);
  iter=0;
  % check that density at end points is low enough
  lddiff=12; % min difference in log-density between mode and end-points
  minld=ld(minf);
  step=1;
  while minld>(modeld-lddiff)
    minf=minf-step*modes;
    minld=ld(minf);
    iter=iter+1;
    step=step*2;
    if iter>100
      error(['lik_negbin -> init_negbin_norm: ' ...
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
      error(['lik_negbin -> init_negbin_norm: ' ...
             'integration interval maximun not found ' ...
             'even after looking hard!'])
    end
  end
  
  
  function integrand = binomial_norm(f)
  % Logit * Gaussian
    integrand = exp(ldconst + yy*log(1./(1.+exp(-f)))+(N-yy)*log(1-1./(1.+exp(-f)))...
                   - 0.5 * (f-myy_i).^2./sigm2_i);
%     integrand = exp(ldconst ...
%                     +yy*log(x)+(N-yy)*log(1-x) ...
%                     -0.5*(f-myy_i).^2./sigm2_i);
    integrand(isnan(integrand))=0;
  end
  
  function log_int = log_binomial_norm(f)
  % log(Binomial * Gaussian)
  % log_binomial_norm is used to avoid underflow when searching
  % integration interval
  
    log_int = ldconst + yy*log(1./(1.+exp(-f)))+(N-yy)*log(1-1./(1.+exp(-f)))...
                   - 0.5 * (f-myy_i).^2./sigm2_i;
%     log_int = ldconst ...
%               -log(1+exp(-yy.*f)) ...
%               -0.5*(f-myy_i).^2./sigm2_i;
  end
  
  function g = log_binomial_norm_g(f)
  % d/df log(Binomial * Gaussian)
  % derivative of log_logit_norm
    g = -(f-myy_i)./sigm2_i - exp(-f).*(N-yy)./((1+exp(-f)).^2.*(1-1./(1+exp(-f)))) ...
        + exp(-f).*yy./(1+exp(-f));
%     g = yy./(exp(f*yy)+1)...
%         + (myy_i - f)./sigm2_i;
  end
  
  function g2 = log_binomial_norm_g2(f)
  % d^2/df^2 log(Binomial * Gaussian)
  % second derivate of log_logit_norm
    g2 = - (1+exp(2.*f)+exp(f).*(2+N*sigm2_i)./((1+exp(f))^2*sigm2_i));
%     a=exp(f*yy);
%     g2 = -a*(yy./(a+1)).^2 ...
%          -1/sigm2_i;
  end
  
end

function p = lik_binomial_invlink(lik, f, z)
%LIK_BINOMIAL_INVLINK  Returns values of inverse link function
%             
%  Description 
%    P = LIK_BINOMIAL_INVLINK(LIK, F) takes a likelihood structure LIK and
%    latent values F and returns the values of inverse link function P.
%    This subfunction is needed when using gp_predprctmu. 
%
%     See also
%     LIK_BINOMIAL_LL, LIK_BINOMIAL_PREDY
  
  p = logitinv(f);
end

function reclik = lik_binomial_recappend(reclik, ri, lik)
%RECAPPEND  Append the parameters to the record
%
%  Description 
%    RECLIK = GPCF_BINOMIAL_RECAPPEND(RECLIK, RI, LIK) takes a
%    likelihood record structure RECLIK, record index RI and
%    likelihood structure LIK with the current MCMC samples of
%    the parameters. Returns RECLIK which contains all the old
%    samples and the current samples from LIK. This subfunction 
%    is needed when using MCMC sampling (gp_mc).
% 
%  See also
%    GP_MC
  
  if nargin == 2
    reclik.type = 'Binomial';

    % Set the function handles
    reclik.fh.pak = @lik_binomial_pak;
    reclik.fh.unpak = @lik_binomial_unpak;
    reclik.fh.ll = @lik_binomial_ll;
    reclik.fh.llg = @lik_binomial_llg;    
    reclik.fh.llg2 = @lik_binomial_llg2;
    reclik.fh.llg3 = @lik_binomial_llg3;
    reclik.fh.tiltedMoments = @lik_binomial_tiltedMoments;
    reclik.fh.invlink = @lik_binomial_invlink;
    reclik.fh.predprcty = @lik_binomial_predprcty;
    reclik.fh.predy = @lik_binomial_predy;
    reclik.fh.recappend = @likelih_binomial_recappend;
    return
  end

end
