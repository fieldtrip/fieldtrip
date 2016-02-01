function [f, energ, diagn] = esls(f, opt, gp, x, y, z, angle_range)
%ESLS  Markov chain update for a distribution with a Gaussian "prior" 
%      factored out
%
%  Description
%    [F, ENERG, DIAG] = ESLS(F, OPT, GP, X, Y) takes the current
%    latent values F, options structure OPT, Gaussian process
%    structure GP, inputs X and outputs Y. Samples new latent
%    values and returns also energies ENERG and diagnostics DIAG.
%
%    A Markov chain update is applied to the D-element array f leaving a
%    "posterior" distribution
%    P(f) \propto N(f;0,Sigma) L(f)
%    invariant. Where N(0,Sigma) is a zero-mean Gaussian
%    distribution with covariance Sigma. Often L is a likelihood
%    function in an inference problem.
%
%  Reference:
%   Elliptical slice sampling
%   Iain Murray, Ryan Prescott Adams and David J.C. MacKay.
%   The Proceedings of the 13th International Conference on Artificial
%   Intelligence and Statistics (AISTATS), JMLR W&CP 9:541-548, 2010.
%
%  See also
%    GP_MC

% Copyright (c) Iain Murray, September 2009
% Tweak to interface and documentation, September 2010
% Ville Tolvanen, October 2011 - Changed inputs and outputs for the function to
% fit in with other GPstuf samplers. Documentation standardized with other
% GPstuff documentation and modified according to input/output changes.

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

if nargin<=1
  % Return only default options
  if nargin==0
    f=elliptical_sls_opt();
  else
    f=elliptical_sls_opt(f);
  end
  return
end

D = numel(f);
if (nargin < 7) || isempty(angle_range)
  angle_range = 0;
end

if ~isfield(gp.lik, 'nondiagW') || ismember(gp.lik.type, {'LGP' 'LGPC'});
  [K, C] = gp_trcov(gp,x);
else
  if ~isfield(gp.lik,'xtime')
    nl=[0 repmat(size(y,1), 1, length(gp.comp_cf))];
  else
    xtime=gp.lik.xtime;
    nl=[0 size(gp.lik.xtime,1) size(y,1)];
  end
  nl=cumsum(nl);
  nlp=length(nl)-1;
  
  C = zeros(nl(end));
  for i1=1:nlp
    if i1==1 && isfield(gp.lik, 'xtime')
      C((1+nl(i1)):nl(i1+1),(1+nl(i1)):nl(i1+1)) = gp_trcov(gp, xtime, gp.comp_cf{i1});
    else
      C((1+nl(i1)):nl(i1+1),(1+nl(i1)):nl(i1+1)) = gp_trcov(gp, x, gp.comp_cf{i1});
    end
  end
end

if isfield(gp,'meanf')
  [H_m,b_m,B_m]=mean_prep(gp,x,[]);
  C = C + H_m'*B_m*H_m;
end
L=chol(C, 'lower');

cur_log_like = gp.lik.fh.ll(gp.lik, y, f, z);
for i1=1:opt.repeat  

  % Set up the ellipse and the slice threshold
  nu = reshape(L*randn(D, 1), size(f));
  hh = log(rand) + cur_log_like;
  % Set up a bracket of angles and pick a first proposal.
  % "phi = (theta'-theta)" is a change in angle.
  if angle_range <= 0
    % Bracket whole ellipse with both edges at first proposed point
    phi = rand*2*pi;
    phi_min = phi - 2*pi;
    phi_max = phi;
  else
    % Randomly center bracket on current point
    phi_min = -angle_range*rand;
    phi_max = phi_min + angle_range;
    phi = rand*(phi_max - phi_min) + phi_min;
  end
  
  % Slice sampling loop
  slrej = 0;
  while true
    % Compute f for proposed angle difference and check if it's on the slice
    f_prop = f*cos(phi) + nu*sin(phi);
    cur_log_like = gp.lik.fh.ll(gp.lik, y, f_prop, z);
    if (cur_log_like > hh)
      % New point is on slice, ** EXIT LOOP **
      break;
    end
    % Shrink slice to rejected point
    if phi > 0
      phi_max = phi;
    elseif phi < 0
      phi_min = phi;
    else
      error('BUG DETECTED: Shrunk to current position and still not acceptable.');
    end
    % Propose new angle difference
    phi = rand*(phi_max - phi_min) + phi_min;
    slrej = slrej + 1;
  end
  f = f_prop;
end
energ = [];
diagn.rej = slrej;
diagn.opt = opt;
end

function opt = elliptical_sls_opt(opt)
%ELLIPTICAL_SLS_OPT  Default options for elliptical slice sampling
%
%  Description
%    OPT = ELLIPTICAL_SLS_OPT
%      return default options
%    OPT = ELLIPTICAL_SLS_OPT(OPT)
%      fill empty options with default values
%
%  The options and defaults are
%      repeat              - nth accepted value

if nargin < 1
  opt=[];
end

if ~isfield(opt,'repeat')
  opt.repeat=40;
end

end
