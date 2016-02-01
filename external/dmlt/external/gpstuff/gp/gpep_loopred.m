function [Eft, Varft, lpyt, Eyt, Varyt] = gpep_loopred(gp, x, y, varargin)
%GPEP_LOOPRED  Leave-one-out predictions with EP approximation
%
%  Description
%    [EFT, VARFT, LPYT, EYT, VARYT] = GPEP_LOOPRED(GP, X, Y,
%    OPTIONS) takes a Gaussian process structure GP together with a
%    matrix X of training inputs and vector Y of training targets,
%    and evaluates the leave-one-out predictive distribution at
%    inputs X and returns means EFT and variances VARFT of latent
%    variables, the logarithm of the predictive densities PYT, and
%    the predictive means EYT and variances VARYT of observations
%    at input locations X.
%
%    OPTIONS is optional parameter-value pair
%      z      - optional observed quantity in triplet (x_i,y_i,z_i)
%               Some likelihoods may use this. For example, in case of 
%               Poisson likelihood we have z_i=E_i, that is, expected value 
%               for ith case. 
%
%    EP leave-one-out is approximated by leaving-out site-term and
%    using cavity distribution as leave-one-out posterior for the
%    ith latent value. 
%
%  References
%    Manfred Opper and Ole Winther (2000). Gaussian Processes for
%    Classification: Mean-Field Algorithms. Neural Computation,
%    12(11):2655-2684.
%
%    Rasmussen, C. E. and Williams, C. K. I. (2006). Gaussian
%    Processes for Machine Learning. The MIT Press.
%
%  See also
%    GP_LOOPRED, GP_PRED
  
% Copyright (c) 2010-2012  Aki Vehtari, Ville Tolvanen
% Copyright (c) 2011-2012  Ville Tolvanen

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPEP_LOOPRED';
  ip.addRequired('gp', @(x) isstruct(x));
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.parse(gp, x, y, varargin{:});
  z=ip.Results.z;

  [tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,muvec_i,sigm2vec_i,lnZ_i] = ...
      gpep_e(gp_pak(gp), gp, x, y, 'z', z);

  Eft=muvec_i;
  Varft=sigm2vec_i;
  lpyt=lnZ_i;
  if nargout > 3
    [tmp, Eyt, Varyt] = gp.lik.fh.predy(gp.lik, muvec_i, sigm2vec_i, [], z);
  end
  
end
