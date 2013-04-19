function [Ef, Varf, lpy, Ey, Vary] = gpmc_pred(gp, x, y, varargin)
%GPMC_PRED  Predictions with Gaussian Process MCMC approximation.
%
%  Description
%    [EF, VARF] = GPMC_PRED(RECGP, X, Y, XT, OPTIONS) takes a
%    Gaussian processes record structure RECGP (returned by gp_mc)
%    together with a matrix XT of input vectors, matrix X of
%    training inputs and vector Y of training targets. Returns
%    matrices EF and VARF that contain mean and variance of the
%    marginal posterior predictive distribution
%  
%        E[f | xt, y] = E[ E[f | x, y, th] ]
%      Var[f | xt, y] = E[ Var[f | x, y, th] ] + Var[ E[f | x, y, th] ]
%   
%    [EF, VARF, LPY] = GP_PRED(RECGP, X, Y, XT, 'yt', YT, OPTIONS) 
%    returns also the log predictive density LPY of the
%    observations YT at input locations XT
%
%        Lpy = log(p(yt | xt, x, y))
%
%    [EF, VARF, LPY, EY, VARY] = GP_PRED(RECGP, X, Y, XT, OPTIONS) 
%    returns also the predictive means EY and variances VARY for test
%    observations at input locations XT
%
%        Ey(:,i) = E[y | xt, x, y]
%      Vary(:,i) = Var[y | xt, x, y]
%
%    where the latent variables and parameters have been
%    marginalized out.
%
%    [EF, VARF, LPY, EY, VARY] = GPMC_PRED(RECGP, X, Y, OPTIONS)
%    evaluates the predictive distribution at training inputs X
%    and logarithm of the predictive density LPY of the training
%    observations Y.
%  
%    OPTIONS is an optional parameter-value pair
%      predcf - index vector telling which covariance functions are 
%               used for prediction. Default is all (1:gpcfn). See
%               additional information below.
%      tstind - a vector/cell array defining, which rows of X belong 
%               to which training block in *IC type sparse models. 
%               Deafult is []. In case of PIC, a cell array
%               containing index vectors specifying the blocking
%               structure for test data. IN FIC and CS+FIC a vector
%               of length n that points out the test inputs that
%               are also in the training set (if none, set TSTIND=[])
%      yt     - optional observed yt in test points (see below)
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case
%               of Poisson likelihood we have z_i=E_i, that is,
%               expected value for ith case.
%      zt     - optional observed quantity in triplet (xt_i,yt_i,zt_i)
%               Some likelihoods may use this. For example, in case
%               of Poisson likelihood we have z_i=E_i, that is, the
%               expected value for the ith case.
%       
%    NOTE! In case of FIC and PIC sparse approximation the
%    prediction for only some PREDCF covariance functions is just
%    an approximation since the covariance functions are coupled
%    in the approximation and are not strictly speaking additive
%    anymore.
%
%    For example, if you use covariance such as K = K1 + K2 your
%    predictions Ef1 = gpmc_pred(GP, X, Y, X, 'predcf', 1) and Ef2 =
%    gpmc_pred(gp, x, y, x, 'predcf', 2) should sum up to Ef =
%    gpmc_pred(gp, x, y, x). That is Ef = Ef1 + Ef2. With FULL model
%    this is true but with FIC and PIC this is true only
%    approximately. That is Ef \approx Ef1 + Ef2.
%
%    With CS+FIC the predictions are exact if the PREDCF
%    covariance functions are all in the FIC part or if they are
%    CS covariances.
%
%     NOTE! When making predictions with a subset of covariance
%     functions with FIC approximation the predictive variance can
%     in some cases be ill-behaved i.e. negative or unrealistically
%     small. This may happen because of the approximative nature of
%     the prediction.
%
%  See also
%    GPMC_PREDS, GP_PRED, GP_SET, GP_MC
%

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari
  
% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
  
  switch nargout
    case 1
      [Efs] = gpmc_preds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
    case 2
      [Efs, Varfs] = gpmc_preds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
      Varf=mean(Varfs,2) + var(Efs,0,2);
    case 3
      [Efs, Varfs, lpys] = gpmc_preds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
      Varf=mean(Varfs,2) + var(Efs,0,2);
      lpy=log(mean(exp(lpys),2));
    case 4
      [Efs, Varfs, lpys, Eys] = gpmc_preds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
      Varf=mean(Varfs,2) + var(Efs,0,2);
      Ey=mean(Eys,2);
      lpy=log(mean(exp(lpys),2));
    case 5
      [Efs, Varfs, lpys, Eys, Varys] = gpmc_preds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
      Varf=mean(Varfs,2) + var(Efs,0,2);
      Ey=mean(Eys,2);
      Vary=mean(Varys,2) + var(Eys,0,2);
      lpy=log(mean(exp(lpys),2));
  end

end
