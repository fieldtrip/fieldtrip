function [Eft, Varft, lpyt, Eyt, Varyt] = gpia_pred(gp_array, x, y, varargin) 
%GPIA_PRED  Prediction with Gaussian Process GP_IA solution.
%
%  Description
%    [EFT, VARFT] = GPIA_PRED(GP_ARRAY, X, Y, XT, OPTIONS) 
%    takes a cell array of GP structures together with matrix X of
%    training inputs and vector Y of training targets, and
%    evaluates the predictive distribution at test inputs XT with
%    parameters marginalized out with IA. Returns a posterior mean
%    EFT and variance VARFT of latent variables.
%
%    [EFT, VARFT, LPYT] = GPIA_PRED(GP, X, Y, XT, 'yt', YT, OPTIONS)
%    returns also logarithm of the predictive density LPYT of the
%    observations YT at test input locations XT. This can be used
%    for example in the cross-validation. Here Y has to be a vector.
% 
%    [EFT, VARFT, LPYT, EYT, VARYT] = GPIA_PRED(GP, X, Y, XT, OPTIONS)
%    returns also the posterior predictive mean EYT and variance VARYT.
%
%    [EFT, VARFT, LPYT, EYT, VARYT, FT, PFT] = ...
%      GPIA_PRED(GP_ARRAY, X, Y, XT, OPTIONS) 
%    returns also the numerical representation of the marginal
%    posterior of latent variables at each XT. FT is a vector of
%    latent values and PFT_i = p(FT_i) is the posterior density for
%    FT_i.
%
%    [EF, VARF, LPY, EY, VARY] = GPIA_PRED(GP, X, Y, OPTIONS)
%    evaluates the predictive distribution at training inputs X
%    and logarithm of the predictive density LPY of the training
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
% Copyright (c) 2012 Aki Vehtari

% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.    

  ip=inputParser;
  ip.FunctionName = 'GPIA_PRED';
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
      if iscell(gp_array)
        gptype=gp_array{1}.type;
      else
        gptype=gp_array.type;
      end
      switch gptype
        case {'FULL' 'VAR' 'DTC' 'SOR'}
          tstind = [];
        case {'FIC' 'CS+FIC'}
          tstind = 1:size(x,1);
        case 'PIC'
          if iscell(gp_array)
            tstind = gp_array{1}.tr_index;
          else
            tstind = gp_array.tr_index;
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
  if ~isempty(yt);options.yt=yt;end
  if ~isempty(z);options.z=z;end
  if ~isempty(zt);options.zt=zt;end
  if ~isempty(predcf);options.predcf=predcf;end
  if ~isempty(tstind);options.tstind=tstind;end
  
  nGP = numel(gp_array);

  nt = size(xt,1);
  Efts=zeros(nt,nGP); Varfts=zeros(nt,nGP);
  Eyts=[]; Varyts=[];
  Eytts=[]; Varytts=[];
  lpyts=zeros(nt,nGP);
  for i1=1:nGP
    Gp=gp_array{i1};
    P_TH(:,i1) = Gp.ia_weight;
    % make predictions with different models in gp_array
    if nargout > 3
      [Efts(:,i1), Varfts(:,i1), lpytt, Eytts, Varytts]=gp_pred(Gp,x,y,xt,options);
    elseif nargout > 2
      [Efts(:,i1), Varfts(:,i1), lpytt]=gp_pred(Gp,x,y,xt,options);
    else
      [Efts(:,i1), Varfts(:,i1)]=gp_pred(Gp,x,y,xt,options);
    end
    if ~isempty(yt)
      lpyts(:,i1) = lpytt;
    end
    if ~isempty(Eytts)
      Eyts(:,i1) = Eytts;
      Varyts(:,i1) = Varytts;
    end
  end

  % compute combined predictions
  Eft = sum(bsxfun(@times,Efts,P_TH), 2);
  Varft = sum(bsxfun(@times, Varfts, P_TH), 2) + sum(bsxfun(@times,bsxfun(@minus,Efts,Eft).^2, P_TH),2);
  if nargout > 2
    if isempty(yt)
      lpyt = [];
    else
      lpyt = log(sum(bsxfun(@times,exp(lpyts),P_TH),2));
    end
    if nargout > 3
      if isempty(Eyts)
        Eyt = [];
        Varyt = [];
      else        
        Eyt = sum(bsxfun(@times,Eyts,P_TH),2);
        Varyt = sum(bsxfun(@times,Varyts,P_TH),2) + sum(bsxfun(@times,bsxfun(@minus,Eyts,Eyt).^2,P_TH),2);
      end
    end
  end
