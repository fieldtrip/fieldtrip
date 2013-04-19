function ekfcv = gp_kfcve(w, gp, x, y, varargin)
% GPEP_KFCVE Evaluate the mean negative log k-fold-cv predictive 
%            density.
%
%  Description
%    EKFCV = GP_KFCVE(W, GP, X, Y, OPTIONS) takes a parameter
%    vector W, Gaussian process structure GP, a matrix X of input
%    vectors and a matrix Y of targets, and evaluates the mean
%    negative log leave-one-out predictive density using 
%    K-fold cross-validation.
%
%    OPTIONS is optional parameter-value pair
%      z          - optional observed quantity in triplet (x_i,y_i,z_i)
%                   Some likelihoods may use this. For example, in
%                   case of Poisson likelihood we have z_i=E_i,
%                   that is, expected value for ith case.
%      k          - number of folds in CV
%      rstream    - number of a random stream to be used for
%                   permuting the data befor division. This way
%                   same permutation can be obtained for different
%                   models. Default is 1. See doc RandStream for
%                   more information.
%      trindex    - k-fold CV training indices. A cell array with k
%                   fields each containing index vector for respective
%                   training set.
%      tstindex   - k-fold CV test indices. A cell array with k
%                   fields each containing index vector for
%                   respective test set.
%
%  See also
%    GP_KFCV, GP_OPTIM
%

% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
ip.FunctionName = 'GPEP_KFCVE';
ip.addRequired('w', @(x) isempty(x) || ...
               isvector(x) && isreal(x) && all(isfinite(x)));
ip.addRequired('gp', @(x) isstruct(x) || iscell(x));
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('k', 10, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
ip.addParamValue('rstream', 1, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
ip.addParamValue('trindex', [], @(x) ~isempty(x) || iscell(x))
ip.addParamValue('tstindex', [], @(x) ~isempty(x) || iscell(x))
ip.parse(w, gp, x, y, varargin{:});
k=ip.Results.k;
rstream=ip.Results.rstream;
trindex=ip.Results.trindex;
tstindex=ip.Results.tstindex;
z=ip.Results.z;

gp=gp_unpak(gp, w);
crit = gp_kfcv(gp, x, y, 'z', z, 'inf_method', 'fixed',  ...
               'display', 'off', 'k', k, 'rstream', rstream, ...
               'trindex', trindex, 'tstindex', tstindex);
ekfcv=-crit.mlpd_cv;
