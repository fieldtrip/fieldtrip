function lik = lik_inputdependentnoise(varargin)
%lik_inputdependentnoise    Create input-dependent noise likelihood structure 
%
%  Description
%    LIK = LIK_INPUTDEPENDENTNOISE('PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates a Gaussian likelihood with input dependent noise structure 
%    in which the named parameters have the specified values. Any unspecified
%    parameters are set to default values.
%
%    LIK = LIK_INPUTDEPENDENTNOISE(LIK,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a likelihood function structure with the named
%    parameters altered with the specified values.
%
%    Parameters for Gaussian likelihood function [default]
%      sigma2       - variance [0.1]
%      sigma2_prior - prior for sigma2 [prior_logunif]
%      n            - number of observations per input (See using average
%                     observations below)
%
%    Note! If the prior is 'prior_fixed' then the parameter in
%    question is considered fixed and it is not handled in
%    optimization, grid integration, MCMC etc.
%
%    The likelihood is defined as follows:
%                    __ n
%      p(y|f1, f2) = || i=1 N(y_i | f1_i, sigma2*exp(f2_i))
%
%      where f1 is the first latent variable defining the mean of the
%      gaussian distribution, f2 is the second latent variable defining
%      the noise structure and sigma2 is coefficient for noise.
%
%  See also
%    GP_SET, LIK_*, PRIOR_*
%

% Copyright (c) 2007-2010 Jarno Vanhatalo & Jouni Hartikainen
% Copyright (c) 2010 Aki Vehtari
% Copyright (c) 2011 Ville Tolvanen

% This software is distributed under the GNU General Public
% License (version 2 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_INPUTDEPENDENTNOISE';
  ip.addOptional('lik', [], @isstruct);
  ip.addParamValue('sigma2',0.1, @(x) isscalar(x) && x>0);
  ip.addParamValue('sigma2_prior',prior_logunif(), @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  lik=ip.Results.lik;
  
  if isempty(lik)
    init=true;
    lik.nondiagW=true;
    lik.type = 'Inputdependentnoise';
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'Inputdependentnoise')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('sigma2',ip.UsingDefaults)
    lik.sigma2 = ip.Results.sigma2;
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
    lik.fh.pak = @lik_inputdependentnoise_pak;
    lik.fh.unpak = @lik_inputdependentnoise_unpak;
    lik.fh.lp = @lik_inputdependentnoise_lp;
    lik.fh.lpg = @lik_inputdependentnoise_lpg;
    lik.fh.ll = @lik_inputdependentnoise_ll;
    lik.fh.llg = @lik_inputdependentnoise_llg;    
    lik.fh.llg2 = @lik_inputdependentnoise_llg2;
    lik.fh.llg3 = @lik_inputdependentnoise_llg3;
    lik.fh.predy = @lik_inputdependentnoise_predy;
    lik.fh.invlink = @lik_inputdependentnoise_invlink;
    lik.fh.predprcty = @lik_inputdependentnoise_predprcty;
    lik.fh.recappend = @lik_inputdependentnoise_recappend;
  end

  function [w,s] = lik_inputdependentnoise_pak(lik)
  %LIK_INPUTDEPENDENTNOISE_PAK  Combine likelihood parameters into one vector.
  %
  %  Description 
  %    W = LIK_INPUTDEPENDENTNOISE_PAK(LIK) takes a likelihood structure LIK and
  %    combines the parameters into a single row vector W. This is a mandatory 
  %    subfunction used for example in energy and gradient computations.
  %     
  %       w = log(lik.sigma2)
  %
  %   See also
  %   LIK_INPUTDEPENDENTNOISE_UNPAK, GP_PAK
    
  w = []; s = {};
  if isfield(lik.p, 'sigma2') && ~isempty(lik.p.sigma2)
    w = [w log(lik.sigma2)];
    s = [s 'log(gaussian.sigma2)'];
    % Hyperparameters of sigma2
    [wh sh] = lik.p.sigma2.fh.pak(lik.p.sigma2);
    w = [w wh];
    s = [s sh];
  end    
  end


  function [lik, w] = lik_inputdependentnoise_unpak(lik, w)
  %LIK_INPUTDEPENDENTNOISE_UNPAK  Extract likelihood parameters from the vector.
  %
  %  Description
  %    [LIK, W] = LIK_INPUTDEPENDENTNOISE_UNPAK(W, LIK) takes a likelihood
  %    structure LIK and extracts the parameters from the vector W
  %    to the LIK structure. This is a mandatory subfunction used for 
  %    example in energy and gradient computations.
  %     
  %   Assignment is inverse of  
  %       w = log(lik.sigma2)
  %
  %   See also
  %   LIK_INPUTDEPENDENTNOISE_PAK, GP_UNPAK

  if ~isempty(lik.p.sigma2)
    lik.sigma2 = exp(w(1));
    w = w(2:end);
    
    % Hyperparameters of sigma2
    [p, w] = lik.p.sigma2.fh.unpak(lik.p.sigma2, w);
    lik.p.sigma2 = p;
  end
  end


  function lp = lik_inputdependentnoise_lp(lik, varargin)
  %LIK_INPUTDEPENDENTNOISE_LP  log(prior) of the likelihood parameters
  %
  %  Description
  %    LP = LIK_INPUTDEPENDENTNOISE_LP(LIK) takes a likelihood structure 
  %    LIK and returns log(p(th)), where th collects the parameters. This
  %    subfunction is needed when there are likelihood parameters.
  %
  %  See also
  %    LIK_INPUTDEPENDENTNOISE_LLG, LIK_INPUTDEPENDENTNOISE_LLG3, LIK_INPUTDEPENDENTNOISE_LLG2, GPLA_E
    

  % If prior for sigma2 parameter, add its contribution
  lp = 0;

  if ~isempty(lik.p.sigma2)
    likp=lik.p;
    lp = likp.sigma2.fh.lp(lik.sigma2, likp.sigma2) + log(lik.sigma2);
  end
    
  end

  
  function lpg = lik_inputdependentnoise_lpg(lik)
  %LIK_INPUTDEPENDENTNOISE_LPG  d log(prior)/dth of the likelihood 
  %                parameters th
  %
  %  Description
  %    E = LIK_INPUTDEPENDENTNOISE_LPG(LIK) takes a likelihood structure 
  %    LIK and returns d log(p(th))/dth, where th collects the parameters.
  %    This subfunction is needed when there are likelihood parameters.
  %
  %  See also
  %    LIK_INPUTDEPENDENTNOISE_LLG, LIK_INPUTDEPENDENTNOISE_LLG3, LIK_INPUTDEPENDENTNOISE_LLG2, GPLA_G
    

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
  
  function ll = lik_inputdependentnoise_ll(lik, y, ff, z)
  %LIK_INPUTDEPENDENTNOISE_LL  Log likelihood
  %
  %  Description
  %    LL = LIK_INPUTDEPENDENTNOISE_LL(LIK, Y, F, Z) takes a likelihood
  %    structure LIK, incedence counts Y, expected counts Z, and
  %    latent values F. Returns the log likelihood, log p(y|f,z).
  %    This subfunction is needed when using Laplace approximation 
  %    or MCMC for inference with non-Gaussian likelihoods. This
  %    subfunction is also used in information criteria (DIC, WAIC)
  %    computations.
  %
  %  See also
  %    LIK_INPUTDEPENDENTNOISE_LLG, LIK_INPUTDEPENDENTNOISE_LLG3, LIK_INPUTDEPENDENTNOISE_LLG2, GPLA_E
    
    
    f=ff(:);
    n=size(y,1);
    f1=f(1:n);
    f2=f((n+1):2*n);
    expf2=exp(f2);
    expf2(isinf(expf2))=realmax;
    sigma2 = lik.sigma2;
    
    ll = sum(-0.5*log(2*pi.*sigma2.*expf2) - 1./(2.*sigma2.*expf2).*(y-f1).^2);
    
  end

  function llg = lik_inputdependentnoise_llg(lik, y, ff, param, z)
  %LIK_INPUTDEPENDENTNOISE_LLG  Gradient of the log likelihood
  %
  %  Description 
  %    LLG = LIK_INPUTDEPENDENTNOISE_LLG(LIK, Y, F, PARAM) takes a likelihood
  %    structure LIK, incedence counts Y, expected counts Z and
  %    latent values F. Returns the gradient of the log likelihood
  %    with respect to PARAM. At the moment PARAM can be 'param' or
  %    'latent'. This subfunction is needed when using Laplace
  %    approximation or MCMC for inference with non-Gaussian likelihoods.
  %
  %  See also
  %    LIK_INPUTDEPENDENTNOISE_LL, LIK_INPUTDEPENDENTNOISE_LLG2, LIK_INPUTDEPENDENTNOISE_LLG3, GPLA_E


    f=ff(:);
    n=size(y,1);
    f1=f(1:n);
    f2=f((n+1):2*n);
    expf2 = exp(f2);
    sigma2 = lik.sigma2;
    
    switch param
      case 'param'
          
        llg=sum(-0.5./sigma2+(y-f1).^2./(2*expf2.*sigma2^2));
        % correction for the log transformation
        llg = llg.*lik.sigma2;
        
      case 'latent'
        
        llg1= (y-f1)./(sigma2*expf2);
        llg2= -0.5+(y-f1).^2./(2.*expf2.*sigma2);
        
%         llg1=(-y+f1)/(sqrt(2*pi).*(sigma2.*expf2).^(3/2)).*exp(-1/(2.*sigma2*expf2).*(y-f1).^2);
%         llg2=-exp(f2-(y-f1).^2/(2.*expf2.*sigma2)).*sigma2/(2.*sqrt(2.*pi).*(expf2.*sigma2).^(3/2))+exp(-f2-(y-f1).^2/(2*expf2*sigma2)).*(y-f1).^2/(2.*sqrt(2.*pi.*expf2).*sigma2.^(3/2));
        
        llg=[llg1; llg2];
    end
  end

  function [llg2] = lik_inputdependentnoise_llg2(lik, y, ff, param, z)
  %function [pi_vec, pi_mat] = lik_inputdependentnoise_llg2(lik, y, ff, param, z)
  %LIK_INPUTDEPENDENTNOISE_LLG2  Second gradients of the log likelihood
  %
  %  Description        
  %    LLG2 = LIK_INPUTDEPENDENTNOISE_LLG2(LIK, Y, F, PARAM) takes a likelihood
  %    structure LIK, incedence counts Y, expected counts Z, and
  %    latent values F. Returns the Hessian of the log likelihood
  %    with respect to PARAM. At the moment PARAM can be only
  %    'latent'. LLG2 is a vector with diagonal elements of the
  %    Hessian matrix (off diagonals are zero). This subfunction
  %    is needed when using Laplace approximation or EP for
  %    inference with non-Gaussian likelihoods.
  %
  %  See also
  %    LIK_INPUTDEPENDENTNOISE_LL, LIK_INPUTDEPENDENTNOISE_LLG, LIK_INPUTDEPENDENTNOISE_LLG3, GPLA_E


    f=ff(:);
    
    n=length(y);
    f1=f(1:n);
    f2=f((n+1):2*n);
    sigma2 = lik.sigma2;
    expf2=exp(f2);

    switch param
      case 'param'
        
      case 'latent'
        
%         llg2 = [2./(sigma2.*expf2); 3/2.*(y-f1).^2./(sigma2.*expf2)];
%         llg2mat = [diag(1./sqrt(sigma2.*expf2)); diag(-(y-f1)./sqrt(sigma2.*expf2))];
        
        llg2_11=-1./(sigma2.*expf2);
        llg2_12=-(y-f1)./(sigma2*expf2);
        llg2_22=-(y-f1).^2./(2.*sigma2.*expf2);
        
        llg2 = [llg2_11 llg2_12; llg2_12 llg2_22];
        
      case 'latent+param'
        
        llg2_1=-(y-f1)./(expf2.*sigma2^2);
        llg2_2=-(y-f1).^2./(2*sigma2.^2.*expf2);
        
        llg2=[llg2_1; llg2_2];
        % correction for the log transformation
        llg2 = llg2.*lik.sigma2;
        
    end
  end    
  
  function llg3 = lik_inputdependentnoise_llg3(lik, y, ff, param, z)
  %LIK_INPUTDEPENDENTNOISE_LLG3  Third gradients of the log likelihood
  %
  %  Description
  %    LLG3 = LIK_INPUTDEPENDENTNOISE_LLG3(LIK, Y, F, PARAM) takes a likelihood
  %    structure LIK, incedence counts Y, expected counts Z and
  %    latent values F and returns the third gradients of the log
  %    likelihood with respect to PARAM. At the moment PARAM can be
  %    only 'latent'. LLG3 is a vector with third gradients. This
  %    subfunction is needed when using Laplace approximation for
  %    inference with non-Gaussian likelihoods.
  %
  %  See also
  %    LIK_INPUTDEPENDENTNOISE_LL, LIK_INPUTDEPENDENTNOISE_LLG, LIK_INPUTDEPENDENTNOISE_LLG2, GPLA_E, GPLA_G

    f=ff(:);
    
    n=size(y,1);
    f1=f(1:n);
    f2=f((n+1):2*n);
    expf2=exp(f2);
    sigma2 = lik.sigma2;
    
    switch param
      case 'param'
        
      case 'latent'
        nl=2;
        llg3=zeros(nl,nl,nl,n);
        
        % y=0:
        % thrid derivative derivative wrt f1 (11)
        %       llg3(1,1,1,:) = 0
        % thrid derivative derivative wrt f2 (11)
        llg3(1,1,2,:) = 1./(expf2.*sigma2);
        
        % thrid derivative derivative wrt f1 (12/21)
        llg3(1,2,1,:) = 1./(expf2.*sigma2);
        llg3(2,1,1,:) = llg3(1,2,1,:);
        % thrid derivative derivative wrt f2 (12/21)
        llg3(1,2,2,:) = (y-f1)./(expf2.*sigma2);
        llg3(2,1,2,:) = llg3(1,2,2,:);
        
        % thrid derivative derivative wrt f1 (22)
        llg3(2,2,1,:) = llg3(1,2,2,:);
        % thrid derivative derivative wrt f1 (22)
        llg3(2,2,2,:) = (y-f1).^2./(2.*expf2.*sigma2);
      
      case 'latent2+param'
        
        llg3_11=1./(expf2.*sigma2.^2);
        llg3_12=(y-f1)./(expf2.*sigma2.^2);
        llg3_22=(y-f1).^2./(2.*expf2.*sigma2.^2);
        
        
        llg3 = [diag(llg3_11) diag(llg3_12); diag(llg3_12) diag(llg3_22)];
        % correction for the log transformation
        llg3 = llg3.*lik.sigma2;
        
    end
  end
 

  function [lpy, Ey, Vary] = lik_inputdependentnoise_predy(lik, Ef, Varf, yt, z)
  %LIK_INPUTDEPENDENTNOISE_PREDY  Returns the predictive mean, variance and density of y
  %
  %  Description         
  %    [LPY] = LIK_INPUTDEPENDENTNOISE_PREDY(LIK, EF, VARF, YT) takes a
  %    likelihood structure LIK, posterior mean EF, posterior
  %    variance VARF of the latent variable and observations YT and 
  %    returns the logarithm of the predictive density PY of YT, that is 
  %        p(yt | th) = \int p(yt | f, th) p(f|y) df.
  %    This subfunction is needed when computing posterior predictive 
  %    distributions for future observations.
  %        
  %    [LPY, EY, VARY] = LIK_INPUTDEPENDENTNOISE_PREDY(LIK, EF, VARF YT)
  %    Returns also the posterior predictive mean EY and variance VARY of 
  %    the observations related to the latent variables. This subfunction 
  %    is needed when computing posterior predictive distributions for 
  %    future observations.
  %
  %  See also
  %    GPLA_PRED, GPEP_PRED, GPMC_PRED


%     ntest=size(yt,1);
    Ef = Ef(:);
    ntest = 0.5*size(Ef,1);
    Ef1=Ef(1:ntest); Ef2=Ef(ntest+1:end);
    if size(Varf,2) == size(Varf,1)
      Varf1=diag(Varf(1:ntest,1:ntest));Varf2=diag(Varf(ntest+1:end,ntest+1:end));
    else
      if size(Varf,2)==1
        Varf=reshape(Varf,ntest,2);
      end
      Varf1=Varf(:,1); Varf2=Varf(:,2);
    end
    sigma2=lik.sigma2;
    
    lpy = zeros(size(yt));
   
    Ey=zeros(size(yt));
    Vary=zeros(size(yt));
    
    for i2=1:ntest
      m1=Ef1(i2); m2=Ef2(i2);
      s1=sqrt(Varf1(i2)); s2=sqrt(Varf2(i2));
      pd=@(f1,f2) norm_pdf(yt(i2), f1, sqrt(sigma2.*exp(f2))).*norm_pdf(f1,Ef1(i2),sqrt(Varf1(i2))).*norm_pdf(f2,Ef2(i2),sqrt(Varf2(i2)));
      lpy(i2) = log(dblquad(pd, m1-6.*s1, m1+6.*s1, m2-6.*s2, m2+6.*s2));
    end
%     sigma2 = lik.sigma2;
%     Ef1=Ef(1:ntest);
%     Ef2=Ef((ntest+1):2*ntest);
%     Ey = Ef1;
%     Vary = Varf + sigma2.*exp(Ef2);
    

  end

  function prctys = lik_inputdependentnoise_predprcty(lik, Ef, Varf, zt, prcty)
  %LIK_BINOMIAL_PREDPRCTY  Returns the percentiles of predictive density of y
  %
  %  Description
  %    PRCTY = LIK_BINOMIAL_PREDPRCTY(LIK, EF, VARF YT, ZT)
  %    Returns percentiles of the predictive density PY of YT, that is
  %    This requires also the succes counts YT, numbers of trials ZT. This
  %    subfunction is needed when using function gp_predprcty.
  %
  %  See also
  %    GP_PREDPCTY
  
  n=size(Ef,1)./2;
  prcty = prcty./100;
  prcty = norminv(prcty, 0, 1);
  prctys = bsxfun(@plus, Ef(1:n), bsxfun(@times, sqrt(Varf(1:n) + lik.sigma2.*exp(Ef(n+1:end))), prcty));
  
  end

  function p = lik_inputdependentnoise_invlink(lik, f, z)
  %LIK_INPUTDEPENDENTNOISE_INVLINK  Returns values of inverse link function
  %             
  %  Description 
  %    P = LIK_INPUTDEPENDENTNOISE_INVLINK(LIK, F) takes a likelihood structure LIK and
  %    latent values F and returns the values of inverse link function P.
  %    This subfunction is needed when using function gp_predprcty.
  %
  %     See also
  %     LIK_INPUTDEPENDENTNOISE_LL, LIK_INPUTDEPENDENTNOISE_PREDY
  
    p = exp(f);
  end
  
  function reclik = lik_inputdependentnoise_recappend(reclik, ri, lik)
  %RECAPPEND  Append the parameters to the record
  %
  %  Description 
  %    RECLIK = GPCF_INPUTDEPENDENTNOISE_RECAPPEND(RECLIK, RI, LIK) takes a
  %    likelihood record structure RECLIK, record index RI and
  %    likelihood structure LIK with the current MCMC samples of
  %    the parameters. Returns RECLIK which contains all the old
  %    samples and the current samples from LIK. This subfunction 
  %    is needed when using MCMC sampling (gp_mc).
  % 
  %  See also
  %    GP_MC

  % Initialize record
    if nargin == 2
      reclik.type = 'Inputdependentnoise';
      reclik.nondiagW=true;

      % Initialize parameter

      % Set the function handles
      reclik.fh.pak = @lik_inputdependentnoise_pak;
      reclik.fh.unpak = @lik_inputdependentnoise_unpak;
      reclik.fh.lp = @lik_inputdependentnoise_lp;
      reclik.fh.lpg = @lik_inputdependentnoise_lpg;
      reclik.fh.ll = @lik_inputdependentnoise_ll;
      reclik.fh.llg = @lik_inputdependentnoise_llg;    
      reclik.fh.llg2 = @lik_inputdependentnoise_llg2;
      reclik.fh.llg3 = @lik_inputdependentnoise_llg3;
      reclik.fh.predy = @lik_inputdependentnoise_predy;
      reclik.fh.predprcty = @lik_inputdependentnoise_predprcty;
      reclik.fh.invlink = @lik_inputdependentnoise_invlink;
      reclik.fh.recappend = @lik_inputdependentnoise_recappend;
      reclik.p=[];
      if ~isempty(ri.p.sigma2)
        reclik.p.sigma2 = ri.p.sigma2;
      end
      return
    end
    
    reclik.sigma2(ri,:)=lik.sigma2;
    if ~isempty(lik.p.sigma2)
      reclik.p.sigma2 = feval(lik.p.sigma2.fh.recappend, reclik.p.sigma2, ri, lik.p.sigma2);
    end
  end
end
