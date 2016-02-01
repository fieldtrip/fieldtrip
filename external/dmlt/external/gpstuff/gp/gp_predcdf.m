function cdf = gp_predcdf(gp, x, y, varargin) 
%GP_PREDCDF  Predictive CDF evaluated at YT
%
%  Description
%    CDF = GP_PREDCDF(GP, X, Y, XT, 'yt', YT, OPTIONS)
%    takes a GP structure together with matrix X of training
%    inputs and vector Y of training targets, and evaluates the
%    cdf of the predictive distribution at test inputs XT, YT. 
%
%    CDF = GP_PREDCDF(GP, X, Y, OPTIONS) evaluates the
%    cdf of the predictive distribution at training inputs X, Y.
%
%    OPTIONS is optional parameter-value pair
%      tstind - a vector defining, which rows of X belong to which 
%               training block in *IC type sparse models. Default is [].
%               See also GP_PRED.
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case.
%      zt     - optional observed quantity in triplet (xt_i,yt_i,zt_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, the expected 
%               value for the ith case. 
%
%  See also
%    GP_PRED, GP_PAK, GP_UNPAK
%

% Copyright (c) 2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GP_PREDCDF';
  ip.addRequired('gp',@(x) isstruct(x) || iscell(x));
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addOptional('xt', [], @(x) isempty(x) || (isreal(x) && all(isfinite(x(:)))))
  ip.addParamValue('yt', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('zt', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('prct', [5 50 95], @(x) isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('nsamp', 5000, @(x) isreal(x) && all(isfinite(x(:))))
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
  z = ip.Results.z;
  zt = ip.Results.zt;
  prct = ip.Results.prct;
  nsamp = ip.Results.nsamp;
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
  if ~isempty(z);options.z=z;end
  if ~isempty(yt);options.yt=yt;end
  if ~isempty(zt);options.zt=zt;end
  if ~isempty(predcf);options.predcf=predcf;end
  if ~isempty(tstind);options.tstind=tstind;end

  [tn, nin] = size(x);
  
  if iscell(gp) || numel(gp.jitterSigma2)>1 || isfield(gp,'latent_method')
    % gp_array
    if iscell(gp)
      nGP = numel(gp);
      for i1=1:nGP
        Gp=gp{i1};
        P_TH(:,i1)=Gp.ia_weight;
        [Ef,Varf]=gp_pred(Gp,x,y,xt,options);
        cdfs(:,i1)=Gp.lik.fh.predcdf(Gp.lik, Ef, Varf, yt);
      end
      cdf=sum(cdfs.*P_TH, 2);
    elseif numel(gp.jitterSigma2)>1
      % MCMC samples
      [Efs,Varfs]=gpmc_preds(gp,x,y,xt,options);
      nmc=size(gp.jitterSigma2,1);
      for i1=1:nmc
        Gp = take_nth(gp,i1);
        cdfs(:,i1)=Gp.lik.fh.predcdf(Gp.lik, Ef, Varf, yt);
      end
      cdf=mean(cdfs, 2);
    else
      [Ef,Varf]=gp_pred(gp,x,y,xt,options);
      cdf=gp.lik.fh.predcdf(gp.lik, Ef, Varf, yt);
    end
    
  end

end
