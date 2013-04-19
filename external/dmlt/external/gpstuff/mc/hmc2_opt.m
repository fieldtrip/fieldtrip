function opt = hmc2_opt(opt)
%HMC2_OPT  Default options for Hybrid Monte Carlo sampling.
%
%  Description
%    OPT = HMC2_OPT
%      return default options
%    OPT = HMC2_OPT(OPT)
%      fill empty options with default values
%
%  The options and defaults are
%    display (0)
%      1 to display the energy values and rejection threshold at
%        each step of the Markov chain
%      2 to display also position vectors at each step
%    checkgrad (0)
%      1 to check the user defined gradient function
%    steps (10)
%      Defines the trajectory length (i.e. the number of leapfrog
%      steps at each iteration)
%    nsamples (1)
%      the number of samples retained from the Markov chain
%    nomit (0)
%      the number of samples omitted from the start of the chain
%    persistence (0)
%      0 for complete replacement of momentum variables
%      1 if momentum persistence is used
%    decay (0.9)
%      defines the decay used when a persistent update of
%      (leap-frog) momentum is used. Bounded to the interval [0, 1.)
%    stepadj (0.1)
%      the step adjustment used in leap-frogs
%    stepsf ([])
%      the step size function
%    window (1)
%      the size of the acceptance window
%
%  See also
%    HMC2

%	Copyright (c) Aki Vehtari (1998)

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


if nargin < 1
  opt=[];
end

if ~isfield(opt,'display')
  opt.display=0;
end
if ~isfield(opt,'checkgrad')
  opt.checkgrad=0;
end
if ~isfield(opt,'steps') | opt.steps < 1
  opt.steps=10;
end
if ~isfield(opt,'nsamples') | opt.nsamples < 1
  opt.nsamples=1;
end
if ~isfield(opt,'nomit') | opt.nomit < 0
  opt.nomit=0;
end
if ~isfield(opt,'persistence')
  opt.persistence=0;
end
if ~isfield(opt,'decay') | opt.decay < 0 | opt.decay > 1
  opt.decay=0.9;
end
if ~isfield(opt,'stepadj')
  opt.stepadj=0.1;
end
if ~isfield(opt,'stepsf')
  opt.stepsf=[];
end
if ~isfield(opt,'window') | opt.window < 0
  opt.window=1;
end
if opt.window > opt.steps
  opt.window=opt.steps;
end
