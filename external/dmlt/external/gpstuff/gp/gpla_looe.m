function eloo = gpla_looe(w, gp, x, y, varargin)
%GPLA_LOOE  Evaluate the mean negative log leave-one-out predictive 
%           density with Laplace approximation.
%
%  Description
%    LOOE = GPLA_LOOE(W, GP, X, Y, OPTIONS) takes a parameter
%    vector W, Gaussian process structure GP, a matrix X of input
%    vectors and a matrix Y of targets, and evaluates the mean
%    negative log leave-one-out predictive density
%      LOOE = - 1/n sum log p(Y_i | X, Y_{\i}, th) 
%    where th represents the parameters (lengthScale,
%    magnSigma2...), X is inputs and Y is observations.
%
%    Laplace leave-one-out is approximated by using analytic
%    leave-one-out equations for the approximated latent
%    distribution.
%
%    OPTIONS is optional parameter-value pair
%      z - optional observed quantity in triplet (x_i,y_i,z_i)
%          Some likelihoods may use this. For example, in case of
%          Poisson likelihood we have z_i=E_i, that is, expected
%          value for ith case.
%
%  See also
%    GP_LOOE, GPLA_LOOPRED, GPLA_E
%

% Copyright (c) 2011 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.


ip=inputParser;
ip.FunctionName = 'GPLA_LOOE';
ip.addRequired('w', @(x) isempty(x) || ...
               isvector(x) && isreal(x) && all(isfinite(x)));
ip.addRequired('gp', @(x) isstruct(x) || iscell(x));
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.parse(w, gp, x, y, varargin{:});
z=ip.Results.z;

gp=gp_unpak(gp, w);
[tmp, tmp, lpyt] = gpla_loopred(gp, x, y, 'z', z);
eloo=-sum(lpyt);
