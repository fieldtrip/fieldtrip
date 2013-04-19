function [Eft, Varft, lpyt, Eyt, Varyt] = gpia_loopred(gp_array, x, y, varargin)
%GPIA_LOOPRED  Leave-one-out predictions with Gaussian Process IA approximation
%
%  Description
%    [EFT, VARFT, LPYT, EYT, VARYT] = GPIA_LOOPRED(RECGP, X, Y)
%    takes a cell array of GP structures together with matrix X of
%    input vectors, matrix X of training inputs and vector Y of
%    training targets, and evaluates the leave-one-out predictive
%    distribution at inputs X and returns means EFT and variances
%    VARFT of latent variables, the logarithm of the predictive
%    densities PYT, and the predictive means EYT and variances
%    VARYT of observations at input locations X.
%
%    OPTIONS is optional parameter-value pair
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case. 
%      is     - defines if importance sampling weighting is used 'on'
%               (default). If set to 'off', integration points and
%               weights for the full posterior are used.
%  
%    Given Gaussian likelihood or non-Gaussian likelihood and
%    latent method Laplace or EP, LOO-posterior given
%    hyperparameters is computed analytically or with analytic
%    approximation and LOO-posterior of the hyperparameters is
%    approximated using importance sampling reweighted integration
%    points from gp_ia. Optionally integration weights for full data
%    posterior of hyperparameters can be used (by setting option
%    'is' to 'off').
%
%  References:
%    Aki Vehtari and Jouko Lampinen (2002). Bayesian model
%    assessment and comparison using cross-validation predictive
%    densities. Neural Computation, 14(10):2439-2468.
%
%  See also
%   GP_LOOPRED, GP_IA, GP_PRED
%

% Copyright (c) 2009 Ville Pietil�inen
% Copyright (c) 2009-2010 Jarno Vanhatalo
% Copyright (c) 2010,2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPIA_LOOPRED';
  ip.addRequired('gp_array',@iscell);
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('is', 'on', @(x) ismember(x,{'on' 'off'}))
  ip.parse(gp_array, x, y, varargin{:});
  z=ip.Results.z;
  is=ip.Results.is;

  if isfield(gp_array{1},'meanf') & ~isempty(gp_array{1}.meanf)
    error('GPIA_LOOPRED: Mean functions not yet supported');
  end

  nGP = numel(gp_array);
  n=size(x,1);

  for i1=1:nGP
    Gp=gp_array{i1};
    P_TH(1,i1) = Gp.ia_weight;
    % compute leave-one-out predictions for each hyperparameter sample
    % if latent method is MCMC, then these samples are for latent values, too
    if nargout <= 3
      [Efts(:,i1), Varfts(:,i1), lpyts(:,i1)] = gp_loopred(Gp, x, y, 'z', z);
    else
      [Efts(:,i1), Varfts(:,i1), lpyts(:,i1), Eyts(:,i1), Varyts(:,i1)] = gp_loopred(Gp, x, y, 'z', z);
    end
  end

  if isequal(is,'off')
    P_TH=repmat(P_TH,n,1);
  else
    % log importance sampling weights
    lw=-lpyts;
    % normalize weights
    for i2=1:n
      % this works even when lw have large magnitudes
      lw(i2,:)=lw(i2,:)-sumlogs(lw(i2,:));
    end
    % importance sampling weights
    w=exp(lw);
    % check the effective sample size
    m_eff=1./sum(w.^2,2);
    if min(m_eff)<nGP/10
      warning(sprintf('For %d folds the effective sample size in IS is less than m/10',sum(m_eff<(nGP/10))))
    end
  
    % reweight ia weights
    P_TH=bsxfun(@times,P_TH,w);
    P_TH=bsxfun(@rdivide,P_TH,sum(P_TH,2));
  end

  % compute combined predictions
  Eft = sum(Efts.*P_TH, 2);
  Varft = sum(Varfts.*P_TH,2) + sum(bsxfun(@minus,Efts,Eft).^2.*P_TH,2);
  if nargout > 2
    lpyt = log(sum(exp(lpyts+log(P_TH)),2));
    if nargout > 3
      Eyt = sum(Eyts.*P_TH,2);
      Varyt = sum(Varyts.*P_TH,2) + sum(bsxfun(@minus,Eyts,Eyt).^2.*P_TH, 2);
    end
  end
