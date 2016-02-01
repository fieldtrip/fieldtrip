function lik = lik_gaussiansmt(varargin)
%LIK_GAUSSIANSMT  Create a Gaussian scale mixture likelihood structure
%                 with priors producing approximation of the Student's t
%
%  Description
%    LIK = LIK_GAUSSIANSMT('ndata',N,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    creates a scale mixture noise covariance function structure
%    (with priors producing approximation of the Student's t) in
%    which the named parameters have the specified values. Any
%    unspecified parameters are set to default values. Obligatory
%    parameter is 'ndata', which tells the number of data points,
%    that is, number of mixture components.
%
%    LIK = LIK_GAUSSIANSMT(LIK,'PARAM1',VALUE1,'PARAM2,VALUE2,...) 
%    modify a covariance function structure with the named
%    parameters altered with the specified values.
% 
%    Parameters for the Gaussian scale mixture approximation of the
%    Student's t
%      sigma2    - Variances of the mixture components.
%                  The default is 1 x ndata vector of 0.1s.
%      U         - Part of the parameter expansion, see below.
%                  The default is 1 x ndata vector of 1s.
%      tau2      - Part of the parameter expansion, see below.
%                  The default is 0.1.
%      alpha     - Part of the parameter expansion, see below.
%                  The default is 0.5.
%      nu        - Degrees of freedom. The default is 4.
%      nu_prior  - Prior for nu. The default is prior_fixed().
%      gibbs     - Whether Gibbs sampling is 'on' (default) or 'off'.
%
%    Parametrisation and non-informative priors for alpha and tau
%    are same as in Gelman et. al. (2004) page 304-305:
%      y-E[y] ~ N(0, alpha^2 * U), 
%      where U = diag(u_1, u_2, ..., u_n)
%          u_i ~ Inv-Chi^2(nu, tau^2)
%
%    The parameters of this likelihood can be inferred only by
%    Gibbs sampling by calling GP_MC.
%
%    If degrees of freedom nu is given a prior (other than
%    prior_fixed), it is sampled using slice sampling within Gibbs
%    sampling with limits [0,128].
%
%  See also
%    GP_SET, PRIOR_*, LIK_*

% Copyright (c) 1998,1999,2010 Aki Vehtari
% Copyright (c) 2007-2010 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'LIK_GAUSSIANSMT';
  ip.addOptional('lik', [], @isstruct);
  ip.addParamValue('ndata',[], @(x) isscalar(x) && x>0 && mod(x,1)==0);
  ip.addParamValue('sigma2',[], @(x) isvector(x) && all(x>0));
  ip.addParamValue('U',[], @isvector);
  ip.addParamValue('tau2',0.1, @isscalar);
  ip.addParamValue('alpha',0.5, @isscalar);
  ip.addParamValue('nu',4, @isscalar);
  ip.addParamValue('nu_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('censored',[], @(x) isstruct);
  ip.addParamValue('gibbs','on', @(x) ismember(x,{'on' 'off'}));
  ip.parse(varargin{:});
  lik=ip.Results.lik;

  if isempty(lik)
    init=true;
    lik.type = 'Gaussian-smt';
  else
    if ~isfield(lik,'type') || ~isequal(lik.type,'Gaussian-smt')
      error('First argument does not seem to be a valid likelihood function structure')
    end
    init=false;
  end
  
  % Initialize parameters
  if init || ~ismember('ndata',ip.UsingDefaults)
    ndata = ip.Results.ndata;
    lik.ndata=ndata;
    lik.r = zeros(ndata,1);
  end
  if isempty(ndata)
    error('NDATA has to be defined')
  end
  if init || ~ismember('sigma2',ip.UsingDefaults)
    sigma2=ip.Results.sigma2;
    if isempty(sigma2)
      lik.sigma2 = repmat(0.1,ndata,1);
    else
      if (size(sigma2,1) == lik.ndata && size(sigma2,2) == 1)
        lik.sigma2 = sigma2;
      else
        error('The size of sigma2 has to be NDATAx1')
      end
    end
    lik.sigma2 = sigma2;
  end
  if init || ~ismember('U',ip.UsingDefaults)
    U=ip.Results.U;
    if isempty(U)
      lik.U = ones(ndata,1);
    else
      if size(U,1) == lik.ndata
        lik.U = U;
      else
        error('the size of U has to be NDATAx1')
      end
    end
  end
  if init || ~ismember('tau2',ip.UsingDefaults)
    lik.tau2=ip.Results.tau2;
  end
  if init || ~ismember('alpha',ip.UsingDefaults)
    lik.alpha=ip.Results.alpha;
  end
  if init || ~ismember('nu',ip.UsingDefaults)
    lik.nu=ip.Results.nu;
  end
  if init || ~ismember('censored',ip.UsingDefaults)
    censored=ip.Results.censored;
    if ~isempty(censored)
      lik.censored = censored{1};
      yy = censored{2};
      if lik.censored(1) >= lik.censored(2)
        error('lik_gaussiansmt -> if censored model is used, the limits must be given in increasing order.')
      end
      
      imis1 = [];
      imis2 = [];
      if lik.censored(1) > -inf
        imis1 = find(yy<=lik.censored(1));
      end            
      if lik.censored(1) < inf
        imis2 = find(yy>=lik.censored(2));
      end                                
      lik.cy = yy([imis1 ; imis2])';
      lik.imis = [imis1 ; imis2];
    end
  end
  % Initialize prior structure
  lik.p=[];
  lik.p.sigma=[];
  if init || ~ismember('nu_prior',ip.UsingDefaults)
    lik.p.nu=ip.Results.nu_prior;
  end
  % using Gibbs or not
  if init || ~ismember('gibbs',ip.UsingDefaults)
    lik.gibbs = ip.Results.gibbs;
  end
  if init
    % Set the function handles to the subfunctions
    lik.fh.pak = @lik_gaussiansmt_pak;
    lik.fh.unpak = @lik_gaussiansmt_unpak;
    lik.fh.lp = @lik_gaussiansmt_lp;
    lik.fh.lpg = @lik_gaussiansmt_lpg;
    lik.fh.cfg = @lik_gaussiansmt_cfg;
    lik.fh.trcov  = @lik_gaussiansmt_trcov;
    lik.fh.trvar  = @lik_gaussiansmt_trvar;
    lik.fh.gibbs = @lik_gaussiansmt_gibbs;
    lik.fh.recappend = @lik_gaussiansmt_recappend;
  end

end

function [w,s] = lik_gaussiansmt_pak(lik)
  w = []; s = {};
end

function [lik, w] = lik_gaussiansmt_unpak(lik, w)

end

function lp =lik_gaussiansmt_lp(lik)
  lp = 0;
end

function lpg  = lik_gaussiansmt_lpg(lik)
  lpg = [];
end

function DKff  = lik_gaussiansmt_cfg(lik, x, x2)
  DKff = [];
end

function C = lik_gaussiansmt_trcov(lik, x)
%LIK_GAUSSIANSMT_TRCOV  Evaluate training covariance matrix
%                    corresponding to Gaussian noise
%  Description
%    C = LIK_GAUSSIANSMT_TRCOV(GP, TX) takes in covariance function
%    of a Gaussian process GP and matrix TX that contains
%    training input vectors. Returns covariance matrix C. Every
%    element ij of C contains covariance between inputs i and j
%    in TX. This subfunction is needed only in Gaussian likelihoods.
%
%  See also
%    LIK_GAUSSIANSMT_COV, LIK_GAUSSIANSMT_TRVAR, GP_COV, GP_TRCOV

  [n, m] =size(x);
  n1=n+1;
  
  if n ~= lik.ndata
    error(['lik_gaussiansmt -> _trvar: The training variance can be evaluated'... 
           '      only for training data.                                 '])
  end
  
  C = sparse(1:n, 1:n, lik.sigma2, n, n);
end

function C = lik_gaussiansmt_trvar(lik, x)
%LIK_GAUSSIANSMT_TRVAR  Evaluate training variance vector
%                    corresponding to Gaussian noise
%
%  Description
%    C = LIK_GAUSSIANSMT_TRVAR(LIK, TX) takes in covariance function 
%    of a Gaussian process LIK and matrix TX that contains
%    training inputs. Returns variance vector C. Every
%    element i of C contains variance of input i in TX. This
%    subfunction is needed  only in Gaussian likelihoods.
%
%
%  See also
%    LIK_GAUSSIANSMT_COV, GP_COV, GP_TRCOV
  
  [n, m] =size(x);
  if n ~= lik.ndata
    error(['lik_gaussiansmt -> _trvar: The training variance can be evaluated'... 
           '      only for training data.                                 '])
  end
  C = lik.sigma2;
  
end

function [lik, y] = lik_gaussiansmt_gibbs(gp, lik, x, y)
%LIK_GAUSSIANSMT_GIBBS  Function for sampling the sigma2's
%
%  Description
%    Perform Gibbs sampling for the scale mixture variances. This
%    function is likelihood specific.

  [n,m] = size(x);
  
  % Draw a sample of the mean of y. Its distribution is
  % f ~ N(K*inv(C)*y, K - K*inv(C)*K')
  switch gp.type
    case 'FULL'
      sampy = gp_rnd(gp, x, y, x);
    case 'FIC'
      sampy = gp_rnd(gp, x, y, x, 'tstind', 1:n);
    case {'PIC' 'PIC_BLOCK'}
      sampy = gp_rnd(gp, x, y, x, 'tstind', gp.tr_index);
  end
  % Calculate the residual
  r = y-sampy;
  
  U = lik.U;
  t2 = lik.tau2;
  alpha = lik.alpha;
  nu = lik.nu;
  rss2=alpha.^2.*U;
  
  % Perform the gibbs sampling (Gelman et.al. (2004) page 304-305)
  % Notice that 'sinvchi2rand' is parameterized as in Gelman et. al.
  U=sinvchi2rand(nu+1, (nu.*t2+(r./alpha).^2)./(nu+1));        
  shape = n*nu./2;                               % These are parameters...
  invscale = nu.*sum(1./U)./2;                   % used in Gelman et al
  t2=gamrnd(shape, 1./invscale);                 % Notice! The matlab parameterization is different
  alpha2=sinvchi2rand(n,mean(r.^2./U));
  rss2=alpha2.*U;
  if ~isempty(lik.p.nu)
    % Sample nu using Gibbs sampling
    pp = lik.p.nu;
    opt=struct('nomit',4,'display',0,'method','doubling', ...
               'wsize',4,'plimit',5,'unimodal',1,'mmlimits',[0; 128]);
    nu=sls(@(nu) (-sum(sinvchi2_lpdf(U,nu,t2))-pp.fh.lp(nu, pp)),nu,opt);
  end
  lik.sigma2 = rss2;
  lik.U = U;
  lik.tau2 = t2;
  lik.alpha = sqrt(alpha2);
  lik.nu = nu;
  lik.r = r;
  if isfield(lik, 'censored')   
    imis1 = [];
    imis2 = [];
    if lik.censored(1) > -inf
      imis1 = find(y<=lik.censored(1));
      y(imis1)=normrtrand(sampy(imis1),alpha2*U(imis1),lik.censored(1));
    end
    
    if lik.censored(1) < inf
      imis2 = find(y>=lik.censored(2));
      y(imis2)=normltrand(sampy(imis2),alpha2*U(imis2),lik.censored(2));
    end
    lik.cy = y([imis1 ; imis2]);
  end
end

function reccf = lik_gaussiansmt_recappend(reccf, ri, lik)
%RECAPPEND  Record append
%
%  Description
%    RECCF = LIK_GAUSSIANSMT_RECAPPEND(RECCF, RI, LIK)
%    takes a likelihood record structure RECCF, record
%    index RI and likelihood structure LIK with the
%    current MCMC samples of the parameters. Returns
%    RECCF which contains all the old samples and the
%    current samples from LIK . This subfunction is 
%    needed when using MCMC sampling (gp_mc).
%
%  See also
%    GP_MC and GP_MC -> RECAPPEND
  
  
  if nargin == 2
    % Initialize the record
    reccf.type = 'Gaussian-smt';
    lik.ndata = [];
    
    % Initialize parameters
    reccf.sigma2 = [];
    
    % Set the function handles
    reccf.fh.pak = @lik_gaussiansmt_pak;
    reccf.fh.unpak = @lik_gaussiansmt_unpak;
    reccf.fh.lp = @lik_gaussiansmt_lp;
    reccf.fh.lpg = @lik_gaussiansmt_lpg;
    reccf.fh.cfg = @lik_gaussiansmt_cfg;
    reccf.fh.cov = @lik_gaussiansmt_cov;
    reccf.fh.trcov  = @lik_gaussiansmt_trcov;
    reccf.fh.trvar  = @lik_gaussiansmt_trvar;
    reccf.fh.gibbs = @lik_gaussiansmt_gibbs;
    reccf.fh.recappend = @lik_gaussiansmt_recappend;
  else  
    % Append to the record
    reccf.ndata = lik.ndata;
    gpp = lik.p;
  
    % record noiseSigma
    reccf.sigma2(ri,:)=lik.sigma2;
    if ~isempty(lik.nu)
      reccf.nu(ri,:)=lik.nu;
      reccf.U(ri,:) = lik.U;
      reccf.tau2(ri,:) = lik.tau2;
      reccf.alpha(ri,:) = lik.alpha;
      reccf.r(ri,:) = lik.r;
    end
    if isfield(lik, 'censored')
      reccf.cy(ri,:) = lik.cy';
    end
  end
end

