function [Eft, Covft, ljpyt, Eyt, Covyt] = gpia_jpred(gp_array, x, y, varargin)
%GPIA_JPRED  Prediction with Gaussian Process GP_IA solution.
%
%  Description
%    [EFT, COVFT] = GPIA_JPRED(GP_ARRAY, X, Y, XT, OPTIONS) 
%    takes a cell array of GP structures together with matrix X of
%    training inputs and vector Y of training targets, and
%    evaluates the predictive distribution at test inputs XT with
%    parameters marginalized out with IA. Returns a posterior mean
%    EFT and covariance COVFT of latent variables.
%
%    [EFT, COVFT, JPYT, EYT, COVYT] = GPIA_JPRED(GP, X, Y, XT, 'yt', YT, ...)
%    returns also logarithm of the predictive joint density PYT of
%    the observations YT at test input locations XT with parameters
%    marginalized out with IA. This can be used for example in the
%    cross-validation. Here Y has to be vector. Returns also
%    posterior predictive mean EYT and covariance COVYT.
%
%    [EFT, COVFT, LJPYT, EYT, COVYT, FT, PFT] = ...
%      GPIA_JPRED(GP_ARRAY, X, Y, XT, OPTIONS) 
%    returns also the numerical representation of the marginal
%    posterior of latent variables at each XT. FT is a vector of
%    latent values and PFT_i = p(FT_i) is the posterior density for
%    FT_i.
%
%    [EF, COVF, LJPY, EY, COVY, F, PF] = GPIA_JPRED(GP, X, Y, OPTIONS)
%    evaluates the predictive distribution at training inputs X
%    and logarithm of the predictive density PY of the training
%    observations Y.
%  
%    OPTIONS is optional parameter-value pair
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
%    an approximation since the covariance functions are coupled in
%    the approximation and are not strictly speaking additive
%    anymore.
%
%    For example, if you use covariance such as K = K1 + K2 your
%    predictions Eft1 = gpia_pred(gp_array, X, Y, X, 'predcf', 1) and
%    Eft2 = gpia_pred(gp_array, x, y, x, 'predcf', 2) should sum up to
%    Eft = gpia_pred(gp_array, x, y, x). That is Eft = Eft1 + Eft2. With
%    FULL model this is true but with FIC and PIC this is true only
%    approximately. That is Eft \approx Eft1 + Eft2.
%
%    With CS+FIC the predictions are exact if the PREDCF covariance
%    functions are all in the FIC part or if they are CS
%    covariances.
%
%    NOTE! When making predictions with a subset of covariance
%    functions with FIC approximation the predictive variance can
%    in some cases be ill-behaved i.e. negative or unrealistically
%    small. This may happen because of the approximative nature of
%    the prediction.
%
%  See also
%    GP_PRED, GP_SET, GP_IA
%
  
% Copyright (c) 2009 Ville Pietil�inen
% Copyright (c) 2009-2010 Jarno Vanhatalo    
% Copyright (c) 2011-2012 Ville Tolvanen
% Copyright (c) 2012 Aki Vehtari

% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.    

  
  ip=inputParser;
  ip.FunctionName = 'GPIA_JPRED';
  ip.addRequired('gp_array', @iscell);
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addOptional('xt', [], @(x) isempty(x) || (isreal(x) && all(isfinite(x(:)))))
  ip.addParamValue('yt', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('zt', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('predcf', [], @(x) isempty(x) || ...
                   isvector(x) && isreal(x) && all(isfinite(x)&x>0))
  ip.addParamValue('tstind', [], @(x) isempty(x) || iscell(x) ||...
                   (isvector(x) && isreal(x) && all(isfinite(x)&x>0)))
  if numel(varargin)==0 || isnumeric(varargin{1})
    % inputParser should handle this, but it doesn't
    ip.parse(gp_array, x, y, varargin{:});
  else
    ip.parse(gp_array, x, y, [], varargin{:});
  end
  xt=ip.Results.xt;
  yt=ip.Results.yt;
  z=ip.Results.z;
  zt=ip.Results.zt;
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
  
  % pass these forward
  options=struct();
  if ~isempty(ip.Results.yt);options.yt=ip.Results.yt;end
  if ~isempty(ip.Results.z);options.z=ip.Results.z;end
  if ~isempty(ip.Results.zt);options.zt=ip.Results.zt;end
  if ~isempty(ip.Results.predcf);options.predcf=ip.Results.predcf;end
  if ~isempty(ip.Results.tstind);options.tstind=ip.Results.tstind;end
  
  if nargout > 2 && isempty(yt)
    pyt = NaN;
  end
  
  nGP = numel(gp_array);
  
  for i=1:nGP
    P_TH(i,:) = gp_array{i}.ia_weight;
  end

  % Make predictions with different models in gp_array
  for j = 1:nGP
    if isempty(yt)
      [Eft_grid(j,:), Covft_grid(:,:,j)]=gp_jpred(gp_array{j},x,y,xt,options);            
    else
      [Eft_grid(j,:), Covft_grid(:,:,j), ljpyt_grid(j), Eyt_grid(j,:), Covyt_grid(:,:,j)]=gp_jpred(gp_array{j},x,y,xt, options);
    end
  end
  
  % Sample from the combined predictive distribution
  N = 100; % Sample size per gp in gp_array
  pft = zeros(N*size(Eft_grid,1), 1);
  ft = [];
  for j = 1 : size(Eft_grid,1)
    ftt = repmat(Eft_grid(j,:),N,1)+(chol(Covft_grid(:,:,j))'*randn(size(Covft_grid(:,:,j),1),N))';
    ft = [ft; ftt];
  end
  
  % Calculate mean and variance of the distributions
  Eft = mean(ft)';
  Covft = cov(ft);
  
  % Calculate jpyt with weight given in P_TH.
  if nargout > 2
    if ~isempty(yt)
      ljpyt = log(sum(exp(ljpyt_grid)'.*P_TH));
    else
      error('yt must be provided to get ljpyt');
    end
  end
  
  if nargout > 3
    Eyt = sum(Eyt_grid.*repmat(P_TH,1,size(Eyt_grid,2)),1);
    Covyt = zeros(size(Covyt_grid));
    for i = 1:size(Covyt_grid,3)
      Covyt(:,:,i) = Covyt_grid(:,:,i).*P_TH(i) + diag((Eyt_grid(i,:) - Eyt)).^2;
    end
    Covyt = sum(Covyt,3);
  end
  
  if nargout > 6
    for j=1:size(Eft_grid,1)
      for i=1:N
        index = N*(j-1)+i;
        pft(index,:) = mnorm_pdf(ft(index,:), Eft_grid(j,:), Covft_grid(:,:,j));
      end
    end
    % Normalize distributions
    pft = bsxfun(@rdivide,pft,sum(pft,1));
  end
  
