function [Ef, Covf, ljpy, Ey, Covy] = gpmc_jpreds(gp, x, y, varargin)
%GPMC_JPREDS  Predictions with Gaussian Process MCMC approximation.
%
%  Description
%    [EFS, COVFS] = GPMC_JPREDS(RECGP, X, Y, XT, OPTIONS) takes a
%    Gaussian processes record structure RECGP (returned by gp_mc)
%    together with a matrix XT of input vectors, matrix X of
%    training inputs and vector Y of training targets. Returns
%    matrices EFS and COVFS that contain means and covariances of the
%    conditional posterior predictive distributions given RECGP.
%    In case of non-Gaussian likelihood   
%  
%        Efs(:,i) = E[ f(xt) | f_i, th_i, x, y ]
%      Covfs(:,i) = Cov[ f(xt) | f_i, th_i, x, y ]
%
%    and in case of Gaussian likelihood (f integrated analytically)
%  
%        Efs(:,i) = E[ f(xt) | th_i, x, y ]
%      Covfs(:,i) = Cov[ f(xt) | th_i, x, y ]
%  
%    The marginal posterior mean and variance can be evaluated from
%    these as follows (See also GPMC_JPRED):
%  
%        E[f | xt, y] = E[ E[f | x, y, th] ]
%                     = mean(Efs,2)
%      Var[f | xt, y] = E[ Cov[f | x, y, th] ] + Var[ E[f | x, y, th] ]
%                     = mean(Covfs,2) + diag(var(Efs,0,2))
%   
%    [EFS, COVFS, LJPYS] = GPMC_JPREDS(RECGP, X, Y, XT, OPTIONS) 
%    returns also the logarithm of the predictive joint density JPYS.
% 
%        Pys(:,i) = p(yt | xt, x, y, th_i)
%
%    [EFS, COVFS, LJPYS, EYS, COVYS] = 
%      GPMC_JPREDS(RECGP, X, Y, XT, 'yt', YT, OPTIONS) 
%    returns also the posterior predictive means and covariances.
%
%      Eys(:,i) = E[y | xt, x, y, th_i]
%      Covys(:,i) = Cov[y | xt, x, y, th_i]
% 
%    where the latent variables have been marginalized out.
% 
%    [EFS, COVFS, LJPYS, EYS, COVYS] = GPMC_JPREDS(RECGP, X, Y, OPTIONS)
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
%    GPMC_JPRED, GP_JPRED, GP_SET, GP_MC

% Copyright (c) 2007-2010 Jarno Vanhatalo
% Copyright (c) 2011-2012 Ville Tolvanen
% Copyright (c) 2012 Aki Vehtari
  
% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  
  ip=inputParser;
  ip.FunctionName = 'GPMC_JPREDS';
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
  
  tn = size(x,1);
  if nargin < 4
    error('Requires at least 4 arguments');
  end

%   if nargout > 2 && ~isempty(yt)
%     error('mc_pred -> If py is wanted you must provide the vector yt as an optional input.')
%   end
  
  nin  = size(x,2);
  nout = 1;
  nmc=size(gp.jitterSigma2,1);
  
  if isfield(gp, 'latentValues') && ~isempty(gp.latentValues)
    % Non-Gaussian likelihood. The latent variables should be used in
    % place of observations
    y = gp.latentValues';
  else 
    y = repmat(y,1,nmc);
  end

  if strcmp(gp.type, 'PIC_BLOCK') || strcmp(gp.type, 'PIC')
    ind = gp.tr_index;           % block indeces for training points
    gp = rmfield(gp,'tr_index');
  end
  
  % loop over all samples
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
    
    if nargout < 3
      [Ef(:,i1), Covf(:,:,i1)] = gp_jpred(Gp, x, y(:,i1), xt, 'predcf', predcf, 'tstind', tstind);
    else 
      if isfield(gp, 'latentValues')
        
        [Ef(:,i1), Covf(:,:,i1), ljpy(i1), Ey(:,i1), Covy(:,:,i1)] = gp_jpred(Gp, x, y(:,i1), xt, 'predcf', predcf, 'tstind', tstind, 'yt', yt);
        
        if any(diag(Covf(:,:,i1))<=0)
          nzero = find(Covf(:,:,i1)<=0);
          if ~isempty(nzero)
            warning('gp_mc: Some of the Covf diagonal elements are less than or equal to zero. Those are set to 1e-12.') 
            for j=1:length(nzero)
              ni = nzero(j);
              Covf(ni,ni,i1) = 1e-12;    % Ensure positiviness, which may be a problem with FIC
            end
          end
        end
        
      else
        if nargout < 3
          [Ef(:,i1), Covf(:,:,i1)] = gp_jpred(Gp, x, y(:,i1), xt, 'predcf', predcf, 'tstind', tstind);
        else
          [Ef(:,i1), Covf(:,:,i1), ljpy(i1), Ey(:,i1), Covy(:,:,i1)] = gp_jpred(Gp, x, y(:,i1), xt, 'predcf', predcf, 'tstind', tstind, 'yt', yt);
        end
      end            
    end
  end    
end