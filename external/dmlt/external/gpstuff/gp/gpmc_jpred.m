function [Ef, Covf, ljpy, Ey, Covy] = gpmc_jpred(gp, x, y, varargin)
%GPMC_JPRED  Predictions with Gaussian Process MCMC approximation.
%
%  Description
%    [EF, COVF] = GPMC_JPRED(RECGP, X, Y, XT, OPTIONS) takes a
%    Gaussian processes record structure RECGP (returned by gp_mc)
%    together with a matrix XT of input vectors, matrix X of
%    training inputs and vector Y of training targets. Returns
%    matrices EFS and COVFS that contain mean and covariances of the
%    marginal posterior predictive distribution
%  
%        E[f | xt, y] = E[ E[f | x, y, th] ]
%      Var[f | xt, y] = E[ Var[f | x, y, th] ] + Var[ E[f | x, y, th] ]
%   
%    [EF, COVF, LJPY] = GP_JPRED(RECGP, X, Y, XT, 'yt', YT, OPTIONS) 
%    returns also logarithm of the predictive joint density JPY of
%    the observations YT at input locations XT
%
%        Py = p(yt | xt, x, y)
%       
%    [EF, COVF, LJPY, EY, COVY] = GP_JPRED(RECGP, X, Y, XT, OPTIONS) 
%    returns also the predictive means and covariances for test
%    observations at input locations XT
%
%        Ey(:,i) = E[y | xt, x, y]
%      Covy(:,i) = Var[y | xt, x, y]
%
%    where the latent variables and parameters have been
%    marginalized out.
%
%    [EF, COVF, LJPY, EY, COVY] = GPMC_JPRED(RECGP, X, Y, OPTIONS)
%    evaluates the predictive distribution at training inputs X
%    and logarithm of the predictive density PY of the training
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
%     NOTE! In case of FIC and PIC sparse approximation the
%     prediction for only some PREDCF covariance functions is just
%     an approximation since the covariance functions are coupled
%     in the approximation and are not strictly speaking additive
%     anymore.
%
%     For example, if you use covariance such as K = K1 + K2 your
%     predictions Ef1 = gpmc_pred(GP, X, Y, X, 'predcf', 1) and Ef2 =
%     gpmc_pred(gp, x, y, x, 'predcf', 2) should sum up to Ef =
%     gpmc_pred(gp, x, y, x). That is Ef = Ef1 + Ef2. With FULL model
%     this is true but with FIC and PIC this is true only
%     approximately. That is Ef \approx Ef1 + Ef2.
%
%     With CS+FIC the predictions are exact if the PREDCF
%     covariance functions are all in the FIC part or if they are
%     CS covariances.
%
%     NOTE! When making predictions with a subset of covariance
%     functions with FIC approximation the predictive variance can
%     in some cases be ill-behaved i.e. negative or unrealistically
%     small. This may happen because of the approximative nature of
%     the prediction.
%
%  See also
%    GPMC_JPREDS, GP_JPRED, GP_SET, GP_MC
%

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2011-2012 Ville Tolvanen
  
% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
  
  switch nargout % ugly...
    case 1
      [Efs] = gpmc_jpreds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
    case 2
      [Efs, Covfs] = gpmc_jpreds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
      Covf = zeros(size(Covfs(:,:,1)));
      n = size(Covfs,3);
      for i = 1:n
        Covf = Covf + Covfs(:,:,i);
      end
      Covf = Covf./n + diag(var(Efs,0,2));
    case 3
      [Efs, Covfs, ljpys] = gpmc_jpreds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
      Covf = zeros(size(Covfs(:,:,1)));
      n = size(Covfs,3);
      for i = 1:n
        Covf = Covf + Covfs(:,:,i);
      end
      Covf = Covf./n + diag(var(Efs,0,2));
      ljpy=log(mean(exp(ljpys),2));
    case 4
      [Efs, Covfs, ljpys, Eys] = gpmc_jpreds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
      Covf = zeros(size(Covfs(:,:,1)));
      n = size(Covfs,3);
      for i = 1:n
        Covf = Covf + Covfs(:,:,i);
      end
      Covf = Covf./n + diag(var(Efs,0,2));      
      ljpy=log(mean(exp(ljpys),2));
      Ey = mean(Eys,2);
    case 5
      [Efs, Covfs, ljpys, Eys, Covys] = gpmc_jpreds(gp, x, y, varargin{:});
      Ef=mean(Efs,2);
      Covf = zeros(size(Covfs(:,:,1)));
      n = size(Covfs,3);
      for i = 1:n
        Covf = Covf + Covfs(:,:,i);
      end
      Covf = Covf./n + diag(var(Efs,0,2));
      Ey=mean(Eys,2);
      Covy = zeros(size(Covys(:,:,1)));
      n = size(Covys,3);
      for i = 1:n
        Covy = Covy + Covys(:,:,i);
      end
      Covy = Covy./n + diag(var(Eys,0,2));
      ljpy=log(mean(exp(ljpys),2));
  end

end
