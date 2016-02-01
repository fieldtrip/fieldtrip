function x_new = sls1mm(f, x_0, opt, gradf, varargin)
%SLS1MM  Markov Chain Monte Carlo sampling using Slice Sampling
%
%  Description
%    SLS1MM is faster streamlined version of SLS for generating samples
%    from one dimensional distribution with known min-max limits.
%
%    SAMPLE = SLS1MM(F, X, OPTIONS) uses slice sampling to generate
%      single value from the *one dimensional* distribution P ~ EXP(-F), 
%      where F is the first argument to SLS1MM. Markov chain starts from
%      point X and the sampling is made using min-max slice sampling.
%      See SLS1MM_OPT for details.
%
%    SAMPLES = SLS1MM(F, X, OPTIONS, [], P1, P2, ...) allows additional
%      arguments to be passed to F(). The fourth argument is ignored,
%      but included for compatibility with HMC2 and the optimisers.
%
%  See SLS1MM_OPT for the optional parameters in the OPTIONS structure.
%  Note, that unlike in SLS, missing fields give an error.
%
%  See also
%    SLS1MM_OPT, SLS

%  Based on "Slice Sampling" by Radford M. Neal in "The Annals of Statistics"
%  2003, Vol. 31, No. 3, 705-767, (c) Institute of Mathematical Statistics, 2003

%       Copyright (c) 2003-2004 Toni Auranen
%       Copyright (c) 2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Set up some variables
maxiter = opt.maxiter;
l = opt.mmlimits(1);
r = opt.mmlimits(2);

% Generate sample
y = -f(x_0,varargin{:}) + log(rand);
x_new = x_0;
for iter=1:maxiter
  x_new = l + (r-l).*rand;
  y_new = -f(x_new,varargin{:});
  if y < y_new
    return;
  end
  if x_new < x_0
    l = x_new;
  else
    r = x_new;
  end
end
if iter+1 > maxiter
  fprintf('Maximum number (%d) of iterations reached during shrinkage.\n',maxiter);
  error('Check function F, decrease the interval ''mmlimits'' or increase the value of ''maxiter''.');
end
