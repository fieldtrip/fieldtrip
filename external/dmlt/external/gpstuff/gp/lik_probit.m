function lik = lik_probit(varargin)
%LIK_PROBIT  Create a Probit likelihood structure 
%
%  Description
%    LIK = LIK_PROBIT creates Probit likelihood for classification
%    problem with class labels {-1,1}.
%  
%    The likelihood is defined as follows:
%                  __ n
%      p(y|f, z) = || i=1 normcdf(y_i * f_i)
%    
%      where f is the latent value vector.
%
%  See also
%    GP_SET, LIK_*
%

% Copyright (c) 2007      Jaakko Riihim�ki
% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_PROBIT';
  ip.addOptional('lik', [], @isstruct);
  ip.parse(varargin{:});
  lik=ip.Results.lik;

  if isempty(lik)
    init=true;
    lik.type = 'Probit';
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'Probit')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end

  if init
    % Set the function handles to the subfunctions
    lik.fh.pak = @lik_probit_pak;
    lik.fh.unpak = @lik_probit_unpak;
    lik.fh.ll = @lik_probit_ll;
    lik.fh.llg = @lik_probit_llg;    
    lik.fh.llg2 = @lik_probit_llg2;
    lik.fh.llg3 = @lik_probit_llg3;
    lik.fh.tiltedMoments = @lik_probit_tiltedMoments;
    lik.fh.predy = @lik_probit_predy;
    lik.fh.invlink = @lik_probit_invlink;
    lik.fh.recappend = @lik_probit_recappend;
  end

end

function [w,s] = lik_probit_pak(lik)
%LIK_PROBIT_PAK  Combine likelihood parameters into one vector.
%
%  Description 
%    W = LIK_PROBIT_PAK(LIK) takes a likelihood structure LIK and
%    returns an empty verctor W. If Probit likelihood had
%    parameters this would combine them into a single row vector
%    W (see e.g. lik_negbin). This is a mandatory subfunction used 
%    for example in energy and gradient computations.
%       
%     See also
%     LIK_NEGBIN_UNPAK, GP_PAK

  w = []; s = {};
end


function [lik, w] = lik_probit_unpak(lik, w)
%LIK_PROBIT_UNPAK  Extract likelihood parameters from the vector.
%
%  Description
%    W = LIK_PROBIT_UNPAK(W, LIK) Doesn't do anything.
% 
%    If Probit likelihood had parameters this would extracts them
%    parameters from the vector W to the LIK structure. This is a
%    mandatory subfunction used for example in energy and gradient 
%    computations.
%       
%  See also
%    LIK_PROBIT_PAK, GP_UNPAK

  lik=lik;
  w=w;
end

function ll = lik_probit_ll(lik, y, f, z)
%LIK_PROBIT_LL  Log likelihood
%
%  Description
%    E = LIK_PROBIT_LL(LIK, Y, F) takes a likelihood structure
%    LIK, class labels Y, and latent values F. Returns the log
%    likelihood, log p(y|f,z). This subfunction is needed when 
%    using Laplace approximation or MCMC for inference with 
%    non-Gaussian likelihoods. This subfunction is also used 
%    in information criteria (DIC, WAIC) computations.
%
%  See also
%    LIK_PROBIT_LLG, LIK_PROBIT_LLG3, LIK_PROBIT_LLG2, GPLA_E

  if ~isempty(find(abs(y)~=1))
    error('lik_probit: The class labels have to be {-1,1}')
  end
  
  p = y.*f;
  ll = log(norm_cdf(p));
  if any(p<-10)
    % Asymptotic expansion of norm_cdf
    i = find(p<-10);
    c = 1 - 1./p(i).^2.*(1-3./p(i).^2.*(1-5./p(i).^2.*(1-7./p(i).^2)));
    ll(i) = -0.5*log(2*pi)-p(i).^2./2-log(-p(i))+log(c); 
  end
  ll = sum(ll);
end


function llg = lik_probit_llg(lik, y, f, param, z)
%LIK_PROBIT_LLG  Gradient of the log likelihood
%
%  Description
%    LLG = LIK_PROBIT_LLG(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, class labels Y, and latent values F. 
%    Returns the gradient of the log likelihood with respect to
%    PARAM. At the moment PARAM can be 'param' or 'latent'.
%    This subfunction is needed when using Laplace approximation 
%    or MCMC for inference with non-Gaussian likelihoods.
%
%   See also
%   LIK_PROBIT_LL, LIK_PROBIT_LLG2, LIK_PROBIT_LLG3, GPLA_E

  if ~isempty(find(abs(y)~=1))
    error('lik_probit: The class labels have to be {-1,1}')
  end
  
  switch param
    case 'latent'
      p=y.*f;
      ncdf=norm_cdf(p);
      if any(p<-10)
        % Asymptotic expansion of norm_cdf
        i = find(p<-10);
        c = 1 - 1./p(i).^2.*(1-3./p(i).^2.*(1-5./p(i).^2.*(1-7./p(i).^2)));
        ncdf(i) = -0.5*log(2*pi)-p(i).^2./2-log(-p(i))+log(c);
      end
      llg = y.*norm_pdf(f)./ncdf;
  end
end


function llg2 = lik_probit_llg2(lik, y, f, param, z)
%LIK_PROBIT_LLG2  Second gradients of the log likelihood
%
%  Description        
%    LLG2 = LIK_PROBIT_LLG2(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, class labels Y, and latent values F. 
%    Returns the Hessian of the log likelihood with respect to
%    PARAM. At the moment PARAM can be only 'latent'. LLG2 is a
%    vector with diagonal elements of the Hessian matrix (off
%    diagonals are zero). This subfunction is needed when using 
%    Laplace approximation or EP for inference with non-Gaussian 
%    likelihoods.
%
%  See also
%    LIK_PROBIT_LL, LIK_PROBIT_LLG, LIK_PROBIT_LLG3, GPLA_E

  
  if ~isempty(find(abs(y)~=1))
    error('lik_probit: The class labels have to be {-1,1}')
  end
  
  switch param
    case 'latent'
      z = y.*f;
      ncdf=norm_cdf(z);
      if any(z<-10)
        % Asymptotic expansion of norm_cdf
        i = find(z<-10);
        c = 1 - 1./z(i).^2.*(1-3./z(i).^2.*(1-5./z(i).^2.*(1-7./z(i).^2)));
        ncdf(i) = -0.5*log(2*pi)-z(i).^2./2-log(-z(i))+log(c);
      end
      z2 = norm_pdf(f)./ncdf;
      llg2 = -z2.^2 - z.*z2;
  end
end

function llg3 = lik_probit_llg3(lik, y, f, param, z)
%LIK_PROBIT_LLG3  Third gradients of the log likelihood
%
%  Description
%    LLG3 = LIK_PROBIT_LLG3(LIK, Y, F, PARAM) takes a likelihood
%    structure LIK, class labels Y, and latent values F and
%    returns the third gradients of the log likelihood with
%    respect to PARAM. At the moment PARAM can be only 'latent'. 
%    LLG3 is a vector with third gradients. This subfunction is 
%    needed when using Laplace approximation for inference with 
%    non-Gaussian likelihoods.
%
%  See also
%    LIK_PROBIT_LL, LIK_PROBIT_LLG, LIK_PROBIT_LLG2, GPLA_E, GPLA_G

  if ~isempty(find(abs(y)~=1))
    error('lik_probit: The class labels have to be {-1,1}')
  end
  
  switch param
    case 'latent'
      z=y.*f;
      ncdf=norm_cdf(z);
      if any(z<-10)
        % Asymptotic expansion of norm_cdf
        i = find(z<-10);
        c = 1 - 1./z(i).^2.*(1-3./z(i).^2.*(1-5./z(i).^2.*(1-7./z(i).^2)));
        ncdf(i) = -0.5*log(2*pi)-z(i).^2./2-log(-z(i))+log(c);
      end
      z2 = norm_pdf(f)./ncdf;
      llg3 = 2.*y.*z2.^3 + 3.*f.*z2.^2 - z2.*(y-y.*f.^2);
  end
end

function [logM_0, m_1, sigm2hati1] = lik_probit_tiltedMoments(lik, y, i1, sigm2_i, myy_i, z)
%LIK_PROBIT_TILTEDMOMENTS  Returns the marginal moments for EP algorithm
%
%  Description
%    [M_0, M_1, M2] = LIK_PROBIT_TILTEDMOMENTS(LIK, Y, I, S2,
%    MYY) takes a likelihood structure LIK, class labels Y, index
%    I and cavity variance S2 and mean MYY. Returns the zeroth
%    moment M_0, mean M_1 and variance M_2 of the posterior
%    marginal (see Rasmussen and Williams (2006): Gaussian
%    processes for Machine Learning, page 55). This subfunction 
%    is needed when using EP for inference with non-Gaussian 
%    likelihoods.
%
%  See also
%    GPEP_E

% don't check this, because this function is called so often by EP
%  if ~isempty(find(abs(y)~=1))
%    error('lik_probit: The class labels have to be {-1,1}')
%  end

  a=realsqrt(1+sigm2_i);
  zi=y(i1).*myy_i./a;
  
  %normc_zi = 0.5.*erfc(-zi./sqrt(2)); % norm_cdf(zi)
  normc_zi = 0.5.*erfc(-zi./1.414213562373095); % norm_cdf(zi)
  %normp_zi = exp(-0.5.*zi.^2-log(2.*pi)./2); %norm_pdf(zi)
  normp_zi = exp(-0.5.*realpow(zi,2)-0.918938533204673); %norm_pdf(zi)
  m_1=myy_i+(y(i1).*sigm2_i.*normp_zi)./(normc_zi.*a); % muhati1
  sigm2hati1=sigm2_i-(sigm2_i.^2.*normp_zi)./((1+sigm2_i).*normc_zi).*(zi+normp_zi./normc_zi); % sigm2hati1
  logM_0 = reallog(normc_zi);
end

function [lpy, Ey, Vary] = lik_probit_predy(lik, Ef, Varf, yt, zt)
%LIK_PROBIT_PREDY    Returns the predictive mean, variance and density of y
%
%  Description       
%    LPY = LIK_PROBIT_PREDY(LIK, EF, VARF, YT)
%    Returns logarithm of the predictive density PY of YT, that is 
%        p(yt | y) = \int p(yt | f) p(f|y) df.
%    This requires also the class labels YT. This subfunction is 
%    needed when computing posterior predictive distributions for 
%    future observations.
%
%    [LPY, EY, VARY] = LIK_PROBIT_PREDY(LIK, EF, VARF) takes a
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

  if nargout > 1
    py1 = norm_cdf(Ef./sqrt(1+Varf));
    Ey = 2*py1 - 1;

    Vary = 1-Ey.^2;
  end

  if ~isempty(find(abs(yt)~=1))
    error('lik_probit: The class labels have to be {-1,1}')
  end
  lpy=[];
  if ~isempty(yt)
    p = Ef.*yt./sqrt(1+Varf);
    lpy = log(norm_cdf(p));    % Probability p(y_new)
    if any(p<-10)
      % Asymptotic expansion of norm_cdf
      i = find(p<-10);
      c = 1 - 1./p(i).^2.*(1-3./p(i).^2.*(1-5./p(i).^2.*(1-7./p(i).^2)));
      lpy(i) = -0.5*log(2*pi)-p(i).^2./2-log(-p(i))+log(c);
    end
  end
end

function p = lik_probit_invlink(lik, f, z)
%LIK_PROBIT_INVLINK  Returns values of inverse link function
%             
%  Description 
%    P = LIK_PROBIT_INVLINK(LIK, F) takes a likelihood structure LIK and
%    latent values F and returns the values of inverse link function P.
%    This subfunction is needed when using function gp_predprctmu.
%
%     See also
%     LIK_PROBIT_LL, LIK_PROBIT_PREDY
  
  p = norm_cdf(f);
end

function reclik = lik_probit_recappend(reclik, ri, lik)
%RECAPPEND  Append the parameters to the record
%
%  Description 
%    RECLIK = GPCF_PROBIT_RECAPPEND(RECLIK, RI, LIK) takes a
%    likelihood record structure RECLIK, record index RI and
%    likelihood structure LIK with the current MCMC samples of
%    the parameters. Returns RECLIK which contains all the old
%    samples and the current samples from LIK. This subfunction
%    is needed when using MCMC sampling (gp_mc).
% 
%  See also:
%  gp_mc

  if nargin == 2
    reclik.type = 'Probit';

    % Set the function handles
    reclik.fh.pak = @lik_probit_pak;
    reclik.fh.unpak = @lik_probit_unpak;
    reclik.fh.ll = @lik_probit_ll;
    reclik.fh.llg = @lik_probit_llg;    
    reclik.fh.llg2 = @lik_probit_llg2;
    reclik.fh.llg3 = @lik_probit_llg3;
    reclik.fh.tiltedMoments = @lik_probit_tiltedMoments;
    reclik.fh.predy = @lik_probit_predy;
    reclik.fh.invlink = @lik_probit_invlink;
    reclik.fh.recappend = @lik_probit_recappend;
  end
end

