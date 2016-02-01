function [Ef, Varf, lpy, Ey, Vary] = gpmc_preds(gp, x, y, varargin)
%GPMC_PREDS  Predictions with Gaussian Process MCMC approximation.
%
%  Description
%    [EFS, VARFS] = GPMC_PREDS(RECGP, X, Y, XT, OPTIONS) takes a
%    Gaussian processes record structure RECGP (returned by gp_mc)
%    together with a matrix XT of input vectors, matrix X of
%    training inputs and vector Y of training targets. Returns
%    matrices EFS and VARFS that contain means and variances of the
%    conditional posterior predictive distributions given RECGP.
%    In case of non-Gaussian likelihood   
%  
%        Efs(:,i) = E[ f(xt) | f_i, th_i, x, y ]
%      Varfs(:,i) = Var[ f(xt) | f_i, th_i, x, y ]
%
%    and in case of Gaussian likelihood (f integrated analytically)
%  
%        Efs(:,i) = E[ f(xt) | th_i, x, y ]
%      Varfs(:,i) = Var[ f(xt) | th_i, x, y ]
%  
%    The marginal posterior mean and variance can be evaluated from
%    these as follows (See also GPMC_PRED and GP_PRED):
%  
%        E[f | xt, y] = E[ E[f | x, y, th] ]
%                     = mean(Efs,2)
%      Var[f | xt, y] = E[ Var[f | x, y, th] ] + Var[ E[f | x, y, th] ]
%                     = mean(Varfs,2) + var(Efs,0,2)
%   
%    [EFS, VARFS, LPYS] = GPMC_PREDS(RECGP, X, Y, XT, 'yt', YT, OPTIONS) 
%    returns also the log predictive density LPYS of the
%    observations YT at input locations XT given RECGP
%
%        Lpys(:,i) = log(p(yt | xt, x, y, th_i))
%
%    [EFS, VARFS, LPYS, EYS, VARYS] = GPMC_PREDS(RECGP, X, Y, XT, OPTIONS) 
%    returns also the predictive means EYS and variances VARYS for test
%    observations at input locations XT given RECGP
%
%        Eys(:,i) = E[y | xt, x, y, th_i]
%      Varys(:,i) = Var[y | xt, x, y, th_i]
%
%    where the latent variables have been marginalized out.
%
%    [EFS, VARFS, LPYS, EYS, VARYS] = GPMC_PREDS(RECGP, X, Y, OPTIONS)
%    evaluates the predictive distribution at training inputs X
%    and logarithm of the predictive density LPYS of the training
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
%     predictions Ef1 = mc_pred(GP, X, Y, X, 'predcf', 1) and Ef2 =
%     mc_pred(gp, x, y, x, 'predcf', 2) should sum up to Ef =
%     mc_pred(gp, x, y, x). That is Ef = Ef1 + Ef2. With FULL model
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
%    MC_PRED, GP_PRED, GP_SET, GP_MC

% Copyright (c) 2007-2010 Jarno Vanhatalo
  
% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
  
  ip=inputParser;
  ip.FunctionName = 'GPMC_PREDS';
  ip.addRequired('gp',@isstruct);
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addOptional('xt', [], @(x) isempty(x) || (isreal(x) && all(isfinite(x(:)))))
  ip.addParamValue('yt', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('zt', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('predcf', [], @(x) isempty(x) || ...
                   isvector(x) && isreal(x) && all(isfinite(x)&x>0))
  ip.addParamValue('tstind', [], @(x) isempty(x) || iscell(x) ||...
                   (isvector(x) && isreal(x) && all(isfinite(x)&x>0)))
  if numel(varargin)==0 || isnumeric(varargin{1})
    % inputParser should handle this, but it doesn't
    ip.parse(gp, x, y, varargin{:});
  else
    ip.parse(gp, x, y, [], varargin{:});
  end
  xt=ip.Results.xt;
  yt=ip.Results.yt;
  zt=ip.Results.zt;
  z=ip.Results.z;
  predcf=ip.Results.predcf;
  tstind=ip.Results.tstind;
  if isempty(xt)
    xt=x;
    if isempty(tstind)
      if iscell(gp)
        gptype=gp{1}.type;
      else
        gptype=gp.type;
      end
      switch gptype
        case {'FULL' 'VAR' 'DTC' 'SOR'}
          tstind = [];
        case {'FIC' 'CS+FIC'}
          tstind = 1:size(x,1);
        case 'PIC'
          if iscell(gp)
            tstind = gp{1}.tr_index;
          else
            tstind = gp.tr_index;
          end
      end
    end
    if isempty(yt)
      yt=y;
    end
    if isempty(zt)
      zt=z;
    end
  end
  options.yt=yt;
  options.zt=zt;
  options.z=z;
  options.predcf=predcf;
  options.tstind=tstind;
  
  tn = size(x,1);

  if nargout > 2 && isempty(yt)
    error('mc_pred -> If lpy is wanted you must provide the vector yt as an optional input.')
  end
  
  nin  = size(x,2);
  nmc=size(gp.jitterSigma2,1);
  
  if isfield(gp, 'latentValues') && ~isempty(gp.latentValues)
    % Non-Gaussian likelihood with MCMC sampling of latents. The latent
    % variables should be used in place of observations
    lv = gp.latentValues';
  end

  if strcmp(gp.type, 'PIC_BLOCK') || strcmp(gp.type, 'PIC')
    ind = gp.tr_index;           % block indeces for training points
    gp = rmfield(gp,'tr_index');
  end
  
  % loop over all samples
  Ey=[]; Vary=[];
  for i1=1:nmc
    Gp = take_nth(gp,i1);
    if isfield(Gp,'latent_method') && isequal(Gp.latent_method,'MCMC')
      Gp = rmfield(Gp,'latent_method');
    end
    
    switch gp.type            
      case 'FULL' 

      case {'FIC' 'CS+FIC'} 
        % Reformat the inducing inputs 
        u = reshape(Gp.X_u,length(Gp.X_u)/nin,nin);
        Gp.X_u = u;

      case {'PIC' 'PIC_BLOCK'}
        % Reformat the inducing inputs 
        u = reshape(Gp.X_u,length(Gp.X_u)/nin,nin);
        Gp.X_u = u;
        Gp.tr_index = ind;
    end
    
    if isfield(gp, 'latentValues') && ~isempty(gp.latentValues)
      % latent values have been sampled with MCMC
      [Ef(:,i1), Varf(:,i1)] = gp_pred(Gp, x, lv(:,i1), xt, options);
      if any(Varf(:,i1) <= 0)
        % Ensure positiveness, which may be a problem with FIC
        Varf(Varf<=0) = 1e-12; 
        warning('gp_mc: Some of the Varf elements are less than or equal to zero. Those are set to 1e-12.') 
      end
      if nargout >= 4
        [lpy(:,i1), Eyt, Varyt] = Gp.lik.fh.predy(Gp.lik, Ef(:,i1), Varf(:,i1), yt, zt);
        if ~isempty(Eyt)
          Ey(:,i1)=Eyt;
          Vary(:,i1)=Varyt;
        end
      elseif nargout == 3
        lpy(:,i1) = Gp.lik.fh.predy(Gp.lik, Ef(:,i1), Varf(:,i1), yt, zt);
      end
    else 
      % Gaussian likelihood or Laplace/EP for latent values
      if nargout <= 2
        [Ef(:,i1), Varf(:,:,i1)] = gp_pred(Gp, x, y, xt, options);
      elseif nargout <=3
        [Ef(:,i1), Varf(:,:,i1), lpy(:,i1)] = gp_pred(Gp, x, y, xt, options);
      else
        [Ef(:,i1), Varf(:,:,i1), lpy(:,i1), Ey(:,i1), Vary(:,i1)] = gp_pred(Gp, x, y, xt, options); 
      end
    end            
  end    
  Varf=squeeze(Varf);
end
