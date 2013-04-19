function lik = lik_negbin(varargin)
%LIK_NEGBIN  Create a Negative-binomial likelihood structure 
%
%  Description
%    LIK = LIK_NEGBIN('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates Negative-binomial likelihood structure in which the
%    named parameters have the specified values. Any unspecified
%    parameters are set to default values.
%  
%    LIK = LIK_NEGBIN(LIK,'PARAM1',VALUE1,'PARAM2,VALUE2,...)
%    modify a likelihood structure with the named parameters
%    altered with the specified values.
%
%    Parameters for Negative-binomial likelihood [default]
%      disper       - dispersion parameter r [10]
%      disper_prior - prior for disper [prior_logunif]
%  
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc.
%
%    The likelihood is defined as follows:
%                  __ n
%      p(y|f, z) = || i=1 [ (r/(r+mu_i))^r * gamma(r+y_i)
%                           / ( gamma(r)*gamma(y_i+1) )
%                           * (mu/(r+mu_i))^y_i ]
%
%    where mu_i = z_i*exp(f_i) and r is the dispersion parameter.
%    z is a vector of expected mean and f the latent value vector
%    whose components are transformed to relative risk
%    exp(f_i). 
%
%    When using the Negbin likelihood you need to give the vector z
%    as an extra parameter to each function that requires also y. 
%    For example, you should call gpla_e as follows: gpla_e(w, gp,
%    x, y, 'z', z)
%
%  See also
%    GP_SET, LIK_*, PRIOR_*
%

% Copyright (c) 2007-2010 Jarno Vanhatalo & Jouni Hartikainen
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  % inputParser checks the arguments and assigns some default values
  ip=inputParser;
  ip.FunctionName = 'LIK_NEGBIN';
  ip.addOptional('lik', [], @isstruct);
  ip.addParamValue('disper',10, @(x) isscalar(x) && x>0);
  ip.addParamValue('disper_prior',prior_logunif(), @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  lik=ip.Results.lik;
  
  if isempty(lik)
    init=true;
    lik.type = 'Negbin';
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'Negbin')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('disper',ip.UsingDefaults)
    lik.disper = ip.Results.disper;
  end
  % Initialize prior structure
  if init
    lik.p=[];
  end
  if init || ~ismember('disper_prior',ip.UsingDefaults)
    lik.p.disper=ip.Results.disper_prior;
  end
  
  if init
    % Set the function handles to the subfunctions
    lik.fh.pak = @lik_negbin_pak;
    lik.fh.unpak = @lik_negbin_unpak;
    lik.fh.lp = @lik_negbin_lp;
    lik.fh.lpg = @lik_negbin_lpg;
    lik.fh.ll = @lik_negbin_ll;
    lik.fh.llg = @lik_negbin_llg;    
    lik.fh.llg2 = @lik_negbin_llg2;
    lik.fh.llg3 = @lik_negbin_llg3;
    lik.fh.tiltedMoments = @lik_negbin_tiltedMoments;
    lik.fh.siteDeriv = @lik_negbin_siteDeriv;
    lik.fh.predy = @lik_negbin_predy;
    lik.fh.predprcty = @lik_negbin_predprcty;
    lik.fh.invlink = @lik_negbin_invlink;
    lik.fh.recappend = @lik_negbin_recappend;
  end

end

function [w,s] = lik_negbin_pak(lik)
%LIK_NEGBIN_PAK  Combine likelihood parameters into one vector.
%
%  Description 
%    W = LIK_NEGBIN_PAK(LIK) takes a likelihood structure LIK and
%    combines the parameters into a single row vector W. This is a 
%    mandatory subfunction used for example in energy and gradient 
%    computations.
%     
%       w = log(lik.disper)
%
%   See also
%   LIK_NEGBIN_UNPAK, GP_PAK
  
  w=[];s={};
  if ~isempty(lik.p.disper)
    w = log(lik.disper);
    s = [s; 'log(negbin.disper)'];
    [wh sh] = lik.p.disper.fh.pak(lik.p.disper);
    w = [w wh];
    s = [s; sh];
  end
end


function [lik, w] = lik_negbin_unpak(lik, w)
%LIK_NEGBIN_UNPAK  Extract likelihood parameters from the vector.
%
%  Description
%    [LIK, W] = LIK_NEGBIN_UNPAK(W, LIK) takes a likelihood
%    structure LIK and extracts the parameters from the vector W
%    to the LIK structure. This is a mandatory subfunction used 
%    for example in energy and gradient computations.
%     
%   Assignment is inverse of  
%       w = log(lik.disper)
%
%   See also
%   LIK_NEGBIN_PAK, GP_UNPAK

  if ~isempty(lik.p.disper)
    lik.disper = exp(w(1));
    w = w(2:end);
    [p, w] = lik.p.disper.fh.unpak(lik.p.disper, w);
    lik.p.disper = p;
  end
end


function lp = lik_negbin_lp(lik, varargin)
%LIK_NEGBIN_LP  log(prior) of the likelihood parameters
%
%  Description
%    LP = LIK_NEGBIN_LP(LIK) takes a likelihood structure LIK and
%    returns log(p(th)), where th collects the parameters. This
%    subfunction is needed if there are likelihood parameters.
%
%  See also
%    LIK_NEGBIN_LLG, LIK_NEGBIN_LLG3, LIK_NEGBIN_LLG2, GPLA_E
  

% If prior for dispersion parameter, add its contribution
  lp=0;
  if ~isempty(lik.p.disper)
    lp = lik.p.disper.fh.lp(lik.disper, lik.p.disper) +log(lik.disper);
  end
  
end


function lpg = lik_negbin_lpg(lik)
%LIK_NEGBIN_LPG  d log(prior)/dth of the likelihood 
%                parameters th
%
%  Description
%    E = LIK_NEGBIN_LPG(LIK) takes a likelihood structure LIK and
%    returns d log(p(th))/dth, where th collects the parameters.
%    This subfunction is needed if there are likelihood parameters.
%
%  See also
%    LIK_NEGBIN_LLG, LIK_NEGBIN_LLG3, LIK_NEGBIN_LLG2, GPLA_G
  
  lpg=[];
  if ~isempty(lik.p.disper)            
    % Evaluate the gprior with respect to disper
    ggs = lik.p.disper.fh.lpg(lik.disper, lik.p.disper);
    lpg = ggs(1).*lik.disper + 1;
    if length(ggs) > 1
      lpg = [lpg ggs(2:end)];
    end
  end
end  

function ll = lik_negbin_ll(lik, y, f, z)
%LIK_NEGBIN_LL  Log likelihood
%
%  Description
%    LL = LIK_NEGBIN_LL(LIK, Y, F, Z) takes a likelihood
%    structure LIK, incedence counts Y, expected counts Z, and
%    latent values F. Returns the log likelihood, log p(y|f,z).
%    This subfunction is needed when using Laplace approximation 
%    or MCMC for inference with non-Gaussian likelihoods. This 
%    subfunction is also used in information criteria (DIC, WAIC) 
%    computations.
%
%  See also
%    LIK_NEGBIN_LLG, LIK_NEGBIN_LLG3, LIK_NEGBIN_LLG2, GPLA_E
  
  if isempty(z)
    error(['lik_negbin -> lik_negbin_ll: missing z!    '... 
           'Negbin likelihood needs the expected number of    '...
           'occurrences as an extra input z. See, for         '...
           'example, lik_negbin and gpla_e.               ']);
  end

  
  r = lik.disper;
  mu = exp(f).*z;
  ll = sum(r.*(log(r) - log(r+mu)) + gammaln(r+y) - gammaln(r) - gammaln(y+1) + y.*(log(mu) - log(r+mu)));
end

function llg = lik_negbin_llg(lik, y, f, param, z)
%LIK_NEGBIN_LLG  Gradient of the log likelihood
%
%  Description 
%    LLG = LIK_NEGBIN_LLG(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, incedence counts Y, expected counts Z and
%    latent values F. Returns the gradient of the log likelihood
%    with respect to PARAM. At the moment PARAM can be 'param' or
%    'latent'. This subfunction is needed when using Laplace 
%    approximation or MCMC for inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_NEGBIN_LL, LIK_NEGBIN_LLG2, LIK_NEGBIN_LLG3, GPLA_E

  if isempty(z)
    error(['lik_negbin -> lik_negbin_llg: missing z!    '... 
           'Negbin likelihood needs the expected number of    '...
           'occurrences as an extra input z. See, for         '...
           'example, lik_negbin and gpla_e.               ']);
  end

  
  mu = exp(f).*z;
  r = lik.disper;
  switch param
    case 'param'      
      % Derivative using the psi function
      llg = sum(1 + log(r./(r+mu)) - (r+y)./(r+mu) + psi(r + y) - psi(r));
      
      % correction for the log transformation
      llg = llg.*lik.disper;
      
% $$$             % Derivative using sum formulation
% $$$             llg = 0;
% $$$             for i1 = 1:length(y)
% $$$                 llg = llg + log(r/(r+mu(i1))) + 1 - (r+y(i1))/(r+mu(i1));
% $$$                 for i2 = 0:y(i1)-1
% $$$                     llg = llg + 1 / (i2 + r);
% $$$                 end
% $$$             end
% $$$             % correction for the log transformation
% $$$             llg = llg.*lik.disper;
    case 'latent'
      llg = y - (r+y).*mu./(r+mu);
  end
end

function llg2 = lik_negbin_llg2(lik, y, f, param, z)
%LIK_NEGBIN_LLG2  Second gradients of the log likelihood
%
%  Description        
%    LLG2 = LIK_NEGBIN_LLG2(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, incedence counts Y, expected counts Z, and
%    latent values F. Returns the Hessian of the log likelihood
%    with respect to PARAM. At the moment PARAM can be only
%    'latent'. LLG2 is a vector with diagonal elements of the
%    Hessian matrix (off diagonals are zero). This subfunction 
%    is needed when using Laplace approximation or EP for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_NEGBIN_LL, LIK_NEGBIN_LLG, LIK_NEGBIN_LLG3, GPLA_E

  if isempty(z)
    error(['lik_negbin -> lik_negbin_llg2: missing z!   '... 
           'Negbin likelihood needs the expected number of    '...
           'occurrences as an extra input z. See, for         '...
           'example, lik_negbin and gpla_e.               ']);
  end

  
  mu = exp(f).*z;
  r = lik.disper;
  switch param
    case 'param'
      
    case 'latent'
      llg2 = - mu.*(r.^2 + y.*r)./(r+mu).^2;
    case 'latent+param'
      llg2 = (y.*mu - mu.^2)./(r+mu).^2;
      
      % correction due to the log transformation
      llg2 = llg2.*lik.disper;
  end
end    

function llg3 = lik_negbin_llg3(lik, y, f, param, z)
%LIK_NEGBIN_LLG3  Third gradients of the log likelihood
%
%  Description
%    LLG3 = LIK_NEGBIN_LLG3(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, incedence counts Y, expected counts Z and
%    latent values F and returns the third gradients of the log
%    likelihood with respect to PARAM. At the moment PARAM can be
%    only 'latent'. LLG3 is a vector with third gradients. This 
%    subfunction is needed when using Laplace approximation for 
%    inference with non-Gaussian likelihoods.
%
%  See also
%    LIK_NEGBIN_LL, LIK_NEGBIN_LLG, LIK_NEGBIN_LLG2, GPLA_E, GPLA_G

  if isempty(z)
    error(['lik_negbin -> lik_negbin_llg3: missing z!   '... 
           'Negbin likelihood needs the expected number of    '...
           'occurrences as an extra input z. See, for         '...
           'example, lik_negbin and gpla_e.               ']);
  end

  
  mu = exp(f).*z;
  r = lik.disper;
  switch param
    case 'param'
      
    case 'latent'
      llg3 = - mu.*(r.^2 + y.*r)./(r + mu).^2 + 2.*mu.^2.*(r.^2 + y.*r)./(r + mu).^3;
    case 'latent2+param'
      llg3 = mu.*(y.*r - 2.*r.*mu - mu.*y)./(r+mu).^3;
      
      % correction due to the log transformation
      llg3 = llg3.*lik.disper;
  end
end

function [logM_0, m_1, sigm2hati1] = lik_negbin_tiltedMoments(lik, y, i1, sigm2_i, myy_i, z)
%LIK_NEGBIN_TILTEDMOMENTS  Returns the marginal moments for EP
%
%  Description
%    [M_0, M_1, M2] = LIK_NEGBIN_TILTEDMOMENTS(LIK, Y, I, S2,
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
  
%  if isempty(z)
%    error(['lik_negbin -> lik_negbin_tiltedMoments: missing z!'... 
%           'Negbin likelihood needs the expected number of            '...
%           'occurrences as an extra input z. See, for                 '...
%           'example, lik_negbin and gpep_e.                       ']);
%  end
  
  yy = y(i1);
  avgE = z(i1);
  r = lik.disper;
  logM_0=zeros(size(yy));
  m_1=zeros(size(yy));
  sigm2hati1=zeros(size(yy));
  
  for i=1:length(i1)
    % get a function handle of an unnormalized tilted distribution
    % (likelihood * cavity = Negative-binomial * Gaussian)
    % and useful integration limits
    [tf,minf,maxf]=init_negbin_norm(yy(i),myy_i(i),sigm2_i(i),avgE(i),r);
    
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
        error('lik_negbin_tilted_moments: sigm2hati1 >= sigm2_i');
      end
    end
    logM_0(i) = log(m_0);
  end
end

function [g_i] = lik_negbin_siteDeriv(lik, y, i1, sigm2_i, myy_i, z)
%LIK_NEGBIN_SITEDERIV  Evaluate the expectation of the gradient
%                      of the log likelihood term with respect
%                      to the likelihood parameters for EP 
%
%  Description [M_0, M_1, M2] =
%    LIK_NEGBIN_SITEDERIV(LIK, Y, I, S2, MYY, Z) takes a
%    likelihood structure LIK, incedence counts Y, expected
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
    error(['lik_negbin -> lik_negbin_siteDeriv: missing z!'... 
           'Negbin likelihood needs the expected number of        '...
           'occurrences as an extra input z. See, for             '...
           'example, lik_negbin and gpla_e.                   ']);
  end

  yy = y(i1);
  avgE = z(i1);
  r = lik.disper;
  
  % get a function handle of an unnormalized tilted distribution 
  % (likelihood * cavity = Negative-binomial * Gaussian)
  % and useful integration limits
  [tf,minf,maxf]=init_negbin_norm(yy,myy_i,sigm2_i,avgE,r);
  % additionally get function handle for the derivative
  td = @deriv;
  
  % Integrate with quadgk
  [m_0, fhncnt] = quadgk(tf, minf, maxf);
  [g_i, fhncnt] = quadgk(@(f) td(f).*tf(f)./m_0, minf, maxf);
  g_i = g_i.*r;

  function g = deriv(f)
    mu = avgE.*exp(f);
    % Derivative using the psi function
    g = 1 + log(r./(r+mu)) - (r+yy)./(r+mu) + psi(r + yy) - psi(r);
    %      % Derivative using the sum formulation
    %      g = 0;
    %      g = g + log(r./(r+mu)) + 1 - (r+yy)./(r+mu);
    %      for i2 = 0:yy-1
    %        g = g + 1 ./ (i2 + r);
    %      end
  end
end

function [lpy, Ey, Vary] = lik_negbin_predy(lik, Ef, Varf, yt, zt)
%LIK_NEGBIN_PREDY  Returns the predictive mean, variance and density of y
%
%  Description  
%    LPY = LIK_NEGBIN_PREDY(LIK, EF, VARF YT, ZT)
%    Returns logarithm of the predictive density PY of YT, that is 
%        p(yt | zt) = \int p(yt | f, zt) p(f|y) df.
%    This requires also the incedence counts YT, expected counts ZT.
%    This subfunction is needed when computing posterior predictive 
%    distributions for future observations.
%
%    [LPY, EY, VARY] = LIK_NEGBIN_PREDY(LIK, EF, VARF) takes a
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
    error(['lik_negbin -> lik_negbin_predy: missing zt!'... 
           'Negbin likelihood needs the expected number of    '...
           'occurrences as an extra input zt. See, for         '...
           'example, lik_negbin and gpla_e.               ']);
  end

  avgE = zt;
  r = lik.disper;
  
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

        F2 = @(f) exp(log(avgE(i1).*exp(f)+((avgE(i1).*exp(f)).^2/r))+norm_lpdf(f,myy_i,sigm_i));
        EVary(i1) = quadgk(F2,minf,maxf);

        F3 = @(f) exp(2*log(avgE(i1))+2*f+norm_lpdf(f,myy_i,sigm_i));
        VarEy(i1) = quadgk(F3,minf,maxf) - Ey(i1).^2;
    end
    Vary = EVary + VarEy;
  end

  % Evaluate the posterior predictive densities of the given observations
  lpy = zeros(length(yt),1);
  for i1=1:length(yt)
    % get a function handle of the likelihood times posterior
    % (likelihood * posterior = Negative-binomial * Gaussian)
    % and useful integration limits
    [pdf,minf,maxf]=init_negbin_norm(...
      yt(i1),Ef(i1),Varf(i1),avgE(i1),r);
    % integrate over the f to get posterior predictive distribution
    lpy(i1) = log(quadgk(pdf, minf, maxf));
  end
end

function prctys = lik_negbin_predprcty(lik, Ef, Varf, zt, prcty)
%LIK_BINOMIAL_PREDPRCTY  Returns the percentiles of predictive density of y
%
%  Description         
%    PRCTY = LIK_BINOMIAL_PREDPRCTY(LIK, EF, VARF YT, ZT)
%    Returns percentiles of the predictive density PY of YT, that is 
%    This requires also the succes counts YT, numbers of trials ZT.
%    This subfunction is needed when using function gp_predprcty.
%
%  See also 
%    GP_PREDPCTY

  if isempty(zt)
    error(['lik_negbin -> lik_negbin_predprcty: missing zt!'... 
           'Negbin likelihood needs the expected number of    '...
           'occurrences as an extra input zt. See, for         '...
           'example, lik_negbin and gpla_e.               ']);
  end
  
  opt=optimset('TolX',.5,'Display','off');
  nt=size(Ef,1);
  prctys = zeros(nt,numel(prcty));
  prcty=prcty/100;
  r = lik.disper;
  for i1=1:nt
    ci = sqrt(Varf(i1));
    for i2=1:numel(prcty)
      minf = nbininv(prcty(i2),r,r./(r+zt(i1).*exp(Ef(i1)-1.96*ci)));
      maxf = nbininv(prcty(i2),r,r./(r+zt(i1).*exp(Ef(i1)+1.96*ci)));
      a=floor(fminbnd(@(a) (quadgk(@(f) nbincdf(a,r,r./(r+zt(i1).*exp(f))).*norm_pdf(f,Ef(i1),ci),Ef(i1)-6*ci,Ef(i1)+6*ci,'AbsTol',1e-4)-prcty(i2)).^2,minf,maxf,opt));
      if quadgk(@(f) nbincdf(a,r,r./(r+zt(i1).*exp(f))).*norm_pdf(f,Ef(i1),ci),Ef(i1)-6*ci,Ef(i1)+6*ci,'AbsTol',1e-4)<prcty(i2)
        a=a+1;
      end
      prctys(i1,i2)=a;
    end
  end
end

function [df,minf,maxf] = init_negbin_norm(yy,myy_i,sigm2_i,avgE,r)
%INIT_NEGBIN_NORM
%
%  Description
%    Return function handle to a function evaluating
%    Negative-Binomial * Gaussian which is used for evaluating
%    (likelihood * cavity) or (likelihood * posterior) Return
%    also useful limits for integration. This is private function
%    for lik_negbin. This subfunction is needed by subfunctions
%    tiltedMoments, siteDeriv and predy.
%  
%  See also
%    LIK_NEGBIN_TILTEDMOMENTS, LIK_NEGBIN_SITEDERIV,
%    LIK_NEGBIN_PREDY
  
% avoid repetitive evaluation of constant part
  ldconst = -gammaln(r)-gammaln(yy+1)+gammaln(r+yy)...
            - log(sigm2_i)/2 - log(2*pi)/2;
  % Create function handle for the function to be integrated
  df = @negbin_norm;
  % use log to avoid underflow, and derivates for faster search
  ld = @log_negbin_norm;
  ldg = @log_negbin_norm_g;
  ldg2 = @log_negbin_norm_g2;

  % Set the limits for integration
  % Negative-binomial likelihood is log-concave so the negbin_norm
  % function is unimodal, which makes things easier
  if yy==0
    % with yy==0, the mode of the likelihood is not defined
    % use the mode of the Gaussian (cavity or posterior) as a first guess
    modef = myy_i;
  else
    % use precision weighted mean of the Gaussian approximation
    % of the Negative-Binomial likelihood and Gaussian
    mu=log(yy/avgE);
    s2=(yy+r)./(yy.*r);
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
  
  function integrand = negbin_norm(f)
  % Negative-binomial * Gaussian
    mu = avgE.*exp(f);
    integrand = exp(ldconst ...
                    +yy.*(log(mu)-log(r+mu))+r.*(log(r)-log(r+mu)) ...
                    -0.5*(f-myy_i).^2./sigm2_i);
  end
  
  function log_int = log_negbin_norm(f)
  % log(Negative-binomial * Gaussian)
  % log_negbin_norm is used to avoid underflow when searching
  % integration interval
    mu = avgE.*exp(f);
    log_int = ldconst...
              +yy.*(log(mu)-log(r+mu))+r.*(log(r)-log(r+mu))...
              -0.5*(f-myy_i).^2./sigm2_i;
  end
  
  function g = log_negbin_norm_g(f)
  % d/df log(Negative-binomial * Gaussian)
  % derivative of log_negbin_norm
    mu = avgE.*exp(f);
    g = -(r.*(mu - yy))./(mu.*(mu + r)).*mu ...
        + (myy_i - f)./sigm2_i;
  end
  
  function g2 = log_negbin_norm_g2(f)
  % d^2/df^2 log(Negative-binomial * Gaussian)
  % second derivate of log_negbin_norm
    mu = avgE.*exp(f);
    g2 = -(r*(r + yy))/(mu + r)^2.*mu ...
         -1/sigm2_i;
  end
  
end

function mu = lik_negbin_invlink(lik, f, z)
%LIK_NEGBIN_INVLINK  Returns values of inverse link function
%             
%  Description 
%    MU = LIK_NEGBIN_INVLINK(LIK, F) takes a likelihood structure LIK and
%    latent values F and returns the values MU of inverse link function.
%    This subfunction is needed when using function gp_predprctmu.
%
%     See also
%     LIK_NEGBIN_LL, LIK_NEGBIN_PREDY
  
  mu = z.*exp(f);
end

function reclik = lik_negbin_recappend(reclik, ri, lik)
%RECAPPEND  Append the parameters to the record
%
%  Description 
%    RECLIK = GPCF_NEGBIN_RECAPPEND(RECLIK, RI, LIK) takes a
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
    reclik.type = 'Negbin';

    % Initialize parameter
    reclik.disper = [];

    % Set the function handles
    reclik.fh.pak = @lik_negbin_pak;
    reclik.fh.unpak = @lik_negbin_unpak;
    reclik.fh.lp = @lik_negbin_lp;
    reclik.fh.lpg = @lik_negbin_lpg;
    reclik.fh.ll = @lik_negbin_ll;
    reclik.fh.llg = @lik_negbin_llg;    
    reclik.fh.llg2 = @lik_negbin_llg2;
    reclik.fh.llg3 = @lik_negbin_llg3;
    reclik.fh.tiltedMoments = @lik_negbin_tiltedMoments;
    reclik.fh.predy = @lik_negbin_predy;
    reclik.fh.predprcty = @lik_negbin_predprcty;
    reclik.fh.invlink = @lik_negbin_invlink;
    reclik.fh.recappend = @lik_negbin_recappend;
    reclik.p=[];
    reclik.p.disper=[];
    if ~isempty(ri.p.disper)
      reclik.p.disper = ri.p.disper;
    end
  else
    
    % Append to the record
    reclik.disper(ri,:)=lik.disper;
    if ~isempty(lik.p)
      reclik.p.disper = lik.p.disper.fh.recappend(reclik.p.disper, ri, lik.p.disper);
    end
  end
end
