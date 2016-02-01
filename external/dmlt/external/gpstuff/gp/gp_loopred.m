function [Eft, Varft, lpyt, Eyt, Varyt] = gp_loopred(gp, x, y, varargin)
%GP_LOOPRED  Leave-one-out predictions assuming Gaussian observation model
%
%  Description
%    [EFT, VARFT, LPYT, EYT, VARYT] = GP_LOOPRED(GP, X, Y, OPTIONS)
%    takes a Gaussian process structure GP, a matrix X of input
%    vectors and a matrix Y of targets, and evaluates the
%    leave-one-out predictive distributions. Returns a posterior
%    mean EFT and variance VARFT of latent variables, the logarithm
%    of the posterior predictive density PYT, and the posterior
%    predictive mean EYT and variance VARYT of observations at
%    input locations X.
%
%    OPTIONS is optional parameter-value pair
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case. 
%      is     - defines if importance sampling weighting is used 'on'
%               (default). If set to 'off', integration points and
%               weights for the full posterior are used. This
%               option is meaningful only for MCMC sample of GPs
%               from GP_MC or array of GPs from GP_IA.
%
%    Given Gaussian likelihood or non-Gaussian likelihood and
%    latent method Laplace or EP, LOO-posterior given
%    hyperparameters is computed analytically or with analytic
%    approximation and possible LOO-posterior of the
%    hyperparameters is approximated using importance sampling
%    (ignored for GP with single value of hyperparameters). 
%    Optionally full data posterior for hyperparameters can be used
%    (by setting option 'is' to 'off'). Given non-Gaussian
%    likelihood and latent method MCMC, LOO-posterior of the
%    hyperparameters and latent values is approximated using
%    importance sampling.
%
%  References:
%    S. Sundararajan and S. S. Keerthi (2001). Predictive
%    Approaches for Choosing Hyperparameters in Gaussian Processes. 
%    Neural Computation 13:1103-1118.
%
%    Aki Vehtari and Jouko Lampinen (2002). Bayesian model
%    assessment and comparison using cross-validation predictive
%    densities. Neural Computation, 14(10):2439-2468.
%
%    Manfred Opper and Ole Winther (2000). Gaussian Processes for
%    Classification: Mean-Field Algorithms. Neural Computation,
%    12(11):2655-2684.
%
%  See also
%   GP_PRED, GPMC_LOOPRED, GPLA_LOOPRED, GPEP_LOOPRED, GPIA_LOOPRED
%

% Copyright (c) 2008-2010 Jarno Vanhatalo
% Copyright (c) 2010,2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

if iscell(gp) || numel(gp.jitterSigma2)>1 || isfield(gp,'latent_method')
  % use inference specific methods
  if iscell(gp)
    fh_loopred=@gpia_loopred;
  elseif numel(gp.jitterSigma2)>1
    fh_loopred=@gpmc_loopred;
  elseif isfield(gp,'latent_method')
    latent_method=gp.latent_method;
    gplik=gp.lik;
    if strcmp(gp.latent_method,'MCMC')
      % single MCMC sample from the posterior does not allow LOO computation
      % LOO using several MCMC samples is done in gpmc_loopred
      fh_loopred=@gp_pred;
    else
      fh_loopred=gp.fh.loopred;
    end
  else
    error('Logical error by coder of this function!')
  end
  switch nargout
    case 1
      [Eft] = fh_loopred(gp, x, y, varargin{:});
    case 2
      [Eft, Varft] = fh_loopred(gp, x, y, varargin{:});
    case 3
      [Eft, Varft, lpyt] = fh_loopred(gp, x, y, varargin{:});
    case 4
      [Eft, Varft, lpyt, Eyt] = fh_loopred(gp, x, y, varargin{:});
    case 5
      [Eft, Varft, lpyt, Eyt, Varyt] = fh_loopred(gp, x, y, varargin{:});
  end
  return
end

% Nothing to parse, but check the arguments anyway
ip=inputParser;
ip.FunctionName = 'GP_LOOPRED';
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.parse(gp, x, y);

if isfield(gp,'meanf') && ~isempty(gp.meanf)
  error('GP_LOOPRED: Mean functions not yet supported');
end

% For single Gaussian process with Gaussian likelihood LOO
% predictions can be computed analytically.
% S. Sundararajan and S. S. Keerthi (2001). Predictive Approaches
% for Choosing Hyperparameters in Gaussian Processes. Neural
% Computation 13:1103-1118.

% gp_looprep returns b=C\y and iCv=diag(inv(C))
% using efficient computation for CS, FIC, PIC, and CS+FIC
[b,iCv]=gp_looprep(gp,x,y);

% LOO-predictions
myy = y - b./iCv;
sigma2 = 1./iCv;

Eft = myy;
Varft = sigma2-gp.lik.sigma2;
Eyt = myy;
Varyt = sigma2;
lpyt = (-0.5 * (log(2*pi) + log(sigma2) + (y-myy).^2./sigma2));

end
