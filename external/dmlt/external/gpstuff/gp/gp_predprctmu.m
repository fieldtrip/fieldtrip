function prctmus = gp_predprctmu(gp, x, y, varargin) 
%GP_PREPRCTMU  Percentiles of the distribution of the location parameter
%
%  Description
%    PRCTMU = GP_PREDPRCTMU(GP, X, Y, XT, OPTIONS)
%    takes a GP structure together with matrix X of training inputs
%    and vector Y of training targets, and evaluates the
%    percentiles of the distribution of the location parameter at
%    test inputs XT.
%
%    PRCTMU = GP_PREDPRCTMU(GP, X, Y, OPTIONS)
%    evaluates the percentiles of the distribution of the location
%    parameter at training inputs X.
%
%    OPTIONS is optional parameter-value pair
%      prct   - percentiles to be computed (default = [5 50 95])
%      nsamp  - determines the number of samples used by GP_RND in case of 
%               MCMC or IA (default = 5000).
%      predcf - index vector telling which covariance functions are 
%               used for prediction. Default is all (1:gpcfn)
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

% Copyright (c) 2011 Ville Tolvanen
% Copyright (c) 2011-2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GP_PREDPRCTMU';
  ip.addRequired('gp',@(x) isstruct(x) || iscell(x));
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addOptional('xt', [], @(x) isempty(x) || (isreal(x) && all(isfinite(x(:)))))
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
    if isempty(zt)
      zt=z;
    end
  end

  % pass these forward
  options=struct();
  if ~isempty(z);options.z=z;end
  if ~isempty(zt);options.zt=zt;end
  if ~isempty(predcf);options.predcf=predcf;end
  if ~isempty(tstind);options.tstind=tstind;end
  
  [tn, nin] = size(x);
  
  if iscell(gp) || numel(gp.jitterSigma2)>1
    % gp_array or MCMC samples
    % the combined latent posterior is not Gaussian
    sampft = gp_rnd(gp, x, y, xt, 'nsamp', nsamp, options);
    if iscell(gp)
      % in the next step we need just one gp from the array
      gp=gp{1};
    end
    if isfield(gp.lik.fh,'trcov')
      % Gaussian likelihood
      prctmus = prctile(sampft, prct, 2);
    else
      prctmus = prctile(gp.lik.fh.invlink(gp.lik, sampft, zt), prct, 2);
    end
  else  
    % single GP 
    % the latent posterior is Gaussian
    [Eft, Varft] = gp_pred(gp, x, y, xt, 'tstind', tstind, options);
    prct = prct./100;
    prct = norminv(prct, 0, 1);
    if isfield(gp.lik.fh,'trcov')
      % Gaussian likelihood
      prctmus = bsxfun(@plus, Eft, bsxfun(@times, sqrt(Varft), prct));
    else
      % Non-Gaussian likelihood
      np = length(prct);
      prctmus = zeros(size(Eft,1),np);
      for i=1:np
        prctmus(:,i) = gp.lik.fh.invlink(gp.lik, Eft+prct(i).*sqrt(Varft), zt);
      end
    end
  end

end
