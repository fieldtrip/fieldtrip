function [Eft1, Eft2, Covf, lpyt] = pred_coxph(gp, x, y, xt, varargin)
% PRED_COXPH Wrapper for returning useful values for coxph likelihood
%
%  Description
%     [EF1, EF2, COVF, LPYT] = PRED_COXPH(GP,X,Y,XT, OPTIONS)
%     Returns predictions for the two latent processes EF1 and EF2,
%     covariance matrix COVF and logarithms of predictive densities at XT.

% Copyright (c) 2012-2013 Ville Tolvanen

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
ip.FunctionName = 'PRED_COXPH';
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('xt',  @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('yt', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('zt', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('predcf', [], @(x) isempty(x) || ...
                 isvector(x) && isreal(x) && all(isfinite(x)&x>0))
ip.addParamValue('tstind', [], @(x) isempty(x) || iscell(x) ||...
                 (isvector(x) && isreal(x) && all(isfinite(x)&x>0)))
ip.parse(gp, x, y, xt, varargin{:});

if ~strcmp(gp.lik.type, 'Coxph')
  error('Likelihood not Coxph')
end
if nargout > 3
  if numel(gp.jitterSigma2)==1
    [Ef, Covf, lpyt] = gp_pred(gp, x, y, xt, varargin{:});
  else
    [Ef, Covf, lpyt] = gpmc_preds(gp, x, y, xt, varargin{:});
  end
else
  if numel(gp.jitterSigma2)==1
    [Ef, Covf] = gp_pred(gp, x, y, xt, varargin{:});
  else
    [Ef, Covf] = gpmc_preds(gp, x, y, xt, varargin{:});    
  end
end
ntime=size(gp.lik.xtime,1);
if isfield(gp.lik, 'stratificationVariables')
  ind_str=gp.lik.stratificationVariables;
  nf1=ntime.*unique([x(:,ind_str); xt(:,ind_str)], 'rows');
else
  nf1=ntime;
end

Eft1 = Ef(1:nf1,:); Ef(1:nf1,:) = []; Eft2 = Ef;

end

