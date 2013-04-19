function [Ef, Varf, lpy, Ey, Vary] = gpmc_loopreds(gp, x, y, varargin)
%GPMC_LOOPREDS  Leave-one-out predictions with Gaussian Process MCMC approximation
%
%  Description
%    [EFT, VARFT, LPYT, EYT, VARYT] = GPMC_LOOPREDS(RECGP, X, Y)
%    takes a Gaussian processes record structure RECGP (returned by
%    gp_mc) together with a matrix X of training inputs and vector
%    Y of training targets, and evaluates the leave-one-out
%    predictive distribution at inputs X and returns means EFT and
%    variances VARFT of latent variables, the logarithm of the
%    predictive densities PYT, and the predictive means EYT and
%    variances VARYT of observations at input locations X.
%
%    OPTIONS is optional parameter-value pair
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case. 
%    Note:
%      - the hyperparameters, hp, are sampled from the full
%        posterior p(hp|x,y) by gp_mc
%      - in case of Gaussian likelihood or non-Gaussian likelihood 
%        and latent method Laplace/EP, the conditonal LOO-CV 
%        distributions p(f_i | x_\i, y_\i, z_\i, hp_s) are computed
%        for each hyperparameter sample hp_s
%      - use gp_loopred to evaluate p(f_i | x_\i, y_\i, z_\i)
%
%  See also
%   GP_LOOPRED
%

% Copyright (c) 2008-2010,2012 Jarno Vanhatalo
% Copyright (c) 2010,2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% Nothing to parse, but check the arguments anyway
ip=inputParser;
ip.FunctionName = 'GPMC_LOOPREDS';
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.parse(gp, x, y, varargin{:});
z=ip.Results.z;

if isfield(gp,'meanf') & ~isempty(gp.meanf)
  error('GPMC_LOOPREDS: Mean functions not yet supported');
end

nmc=size(gp.jitterSigma2,1);

if isfield(gp, 'latentValues') && ~isempty(gp.latentValues)
  % Non-Gaussian likelihood. The latent variables should be used in
  % place of observations
  lv = gp.latentValues';
end

for i1=1:nmc
  Gp = take_nth(gp,i1);
  if isfield(Gp,'latent_method') && isequal(Gp.latent_method,'MCMC')
    Gp = rmfield(Gp,'latent_method');
  end
  
  if isfield(gp, 'latentValues') && ~isempty(gp.latentValues)
    % latent values have been sampled with MCMC
    [Ef(:,i1), Varf(:,i1)] = gp_loopred(Gp, x, lv(:,i1), varargin{:});
    if nargout >=4
      [lpy(:,i1), Ey(:,i1), Vary(:,i1)] = Gp.lik.fh.predy(Gp.lik, Ef(:,i1), Varf(:,i1), y, z);
    else
      lpy(:,i1) = Gp.lik.fh.predy(Gp.lik, Ef(:,i1), Varf(:,i1), y, z);
    end
  else
    % Gaussian likelihood or Laplace/EP for latent values
    if nargout <= 2
      [Ef(:,i1), Varf(:,i1)] = gp_loopred(Gp, x, y, varargin{:});
    elseif nargout <=3
      [Ef(:,i1), Varf(:,i1), lpy(:,1)] = gp_loopred(Gp, x, y, varargin{:});
    else
      [Ef(:,i1), Varf(:,i1), lpy(:,1), Ey(:,i1), Vary(:,i1)] = gp_loopred(Gp, x, y, varargin{:});
    end
  end
end
