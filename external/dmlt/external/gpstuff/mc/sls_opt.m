function opt = sls_opt(opt)
%SLS_OPT Default options for Slice Sampling
%
%  Default options for SLS
%
%    opt = sls_opt;
%      return default options
%    opt = sls_opt(slsopt);
%      fill empty options with default values
%
%  The options and defaults are (default value)
%    nsamples (1)
%      the number of samples retained from the Markov chain
%    nomit (0)
%      the number of samples omitted from the start of the chain
%    display (1)
%      (0) do not display any information while the chain advances
%      (1) verbose-mode (display only the overflows)
%      (2) display all information while the chain advances (energies,
%          rejections in case overrelaxation, overflows when
%          stepping-out or doubling etc.)
%    method ('stepping')
%      whether to use stepping-out ('stepping'), doubling ('doubling'),
%      ('minmax') to grow the slice or ('multi') and ('multimm') for
%      hyperrectangle multivariate sampling or ('shrnk') and ('covmatch')
%      for covariance-adaptive sampling
%    sigma (1)
%      initial crumb standard deviation for adaptive methods ('shrnk') and
%      ('covmatch')
%    overrelaxation (0)
%      whether to use (1) or not (0) overrelaxed slice sampling
%      (stepping-out with bisection is used)
%    alimit (4)
%      integer limiting the endpoint accuracy to w/(2^a) in bisection
%      (overrelaxation)
%    wsize (1)
%      estimate of the typical slice size w (stepping-out, doubling, shrinkage,
%      multi)
%    mlimit (4)
%      integer for limiting the size of a slice to w*m (stepping-out)
%    maxiter (50)
%      maximum number of iterations for the shrinkage procedure, in case
%      this is exceeded the other parameters might not be properly assigned
%    plimit (2)
%      integer limiting the size of a slice to w*2^p (doubling)
%    unimodal (0)
%      whether the distribution is known to be unimodal (1) or not (0)
%      Note! if the distribution is multimodal, the results are not
%            correct if set to (1)
%    mmlimits ([w-(w*m); w+(w*m)])
%      absolute limits for the slice (minmax, multi, multimm)
%
%  If any of the parameters (alimit, wsize, mlimit, plimit, mmlimits) is a scalar
%  and the parameter is vector format, the parameters are made as vectors. Each
%  of the value can also have its own value when the abovementioned parameters
%  must be given as a vector. Overrelaxation can also be a vector so that
%  overrelaxation is used separately for each variable of a multivariate case.
%
%  See also
%    SLS

%  HUT/LCE, 12/2003
%  (c) Toni Auranen
%
% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

%  Version 1.0, 16/1/2004, TA
%  Version 1.03, 9/3/2004, TA
%   - changed some display-settings
%   - changed the variable initialization
%  Version 1.04, 24/3/2004, TA
%   - overrelaxation separately for each variable in multivariate case
%   - added max_iter parameter for shrinkage
%  Version 1.05, 7/4/2004, TA
%   - added nomit-option
%  Version 1.06, 22/4/2004, TA
%   - added unimodality shortcut for stepping-out and doubling
%  Version 1.06b, 27/4/2004, TA
%   - fixed some bugs
%  Version 1.7, 15/3/2005, TA
%   - added the hyperrectangle multivariate sampling

if nargin < 1
  opt=[];
end

if ~isfield(opt,'nsamples')
  opt.nsamples = 1;
end
if ~isfield(opt,'nomit')
  opt.nomit = 0;
end
if ~isfield(opt,'display')
  opt.display = 1;
end
if ~isfield(opt,'method')
  opt.method = 'stepping';
end
if ~isfield(opt,'overrelaxation')
  opt.overrelaxation = 0;
elseif opt.overrelaxation == 1 && (strcmp(opt.method,'doubling') || strcmp(opt.method,'minmax'))
  opt.method = 'stepping';
end
if ~isfield(opt,'alimit')
  opt.alimit = 4;
end
if ~isfield(opt,'wsize')
  opt.wsize = 1;
end
if ~isfield(opt,'mlimit')
  opt.mlimit = 4;
end
if ~isfield(opt,'maxiter')
  opt.maxiter = 50;
end
if ~isfield(opt,'plimit')
  opt.plimit = 2;
end
if ~isfield(opt,'unimodal')
  opt.unimodal = 0;
end
if ~isfield(opt,'mmlimits')
  opt.mmlimits = [opt.wsize-(opt.wsize*opt.mlimit); opt.wsize+(opt.wsize*opt.mlimit)];
end
if ~isfield(opt, 'sigma') || isempty(opt.sigma)
  opt.sigma = 1;
end
