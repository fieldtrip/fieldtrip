function [x, fs, vs, lambda] = scges(f, x, opt, gradf, varargin)
%SCGES	Scaled conjugate gradient optimization with early stopping.
%
%  Description
%    X = SCGES(F, X, OPT, GRADF, P1,P2,...,P1',P2',...) uses a
%    scaled conjugate gradients algorithm to find a local minimum
%    of the function F(X,P1,P2,...) whose gradient is given by
%    GRADF(X,P1,P2,...). Search is early stopped if value of
%    F(X,P1',P2',...) does not decrease.
%
%    Here X is a row vector and F returns a scalar value. The point
%    at which F has a local minimum is returned as X.
%
%    See SCGES_OPT for the optional parameters in the OPT
%    structure.
%
%    Reference: Vehtari et al (2000). On MCMC sampling in Bayesian
%    MLP neural networks, In Shun-Ichi Amari, C. Lee Giles, Marco
%    Gori, and Vincenzo Piuri, editors, IJCNN'2000: Proceedings of
%    the 2000 International Joint Conference on Neural Networks,
%    volume I, pp. 317-322. IEEE.
%
%  See also SCGES_OPT

% Copyright (c) 1996,1997 Christopher M Bishop, Ian T Nabney
% Copyright (c) 2005 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Set empty options to default values
opt=scges_opt(opt);

% Reference to structures is much slower, so...
niters = opt.maxiter;

display = opt.display;

nparams = length(x);

narg=length(varargin);
varargin1=varargin(1:narg/2);
varargin2=varargin((narg/2+1):end);

%  Check gradients
if opt.checkgrad
  feval('gradcheck', x, f, gradf, varargin1{:});
end

sigma0 = 1.0e-4;
fold = feval(f, x, varargin1{:});	% Initial function value.
vold = feval(f, x, varargin2{:});	% Initial function value.
iter = 0;
gradnew = feval(gradf, x, varargin1{:});% Initial gradient.
gradold = gradnew;
d = - gradnew;				% Initial search direction.
success = 1;                            % Force calculation of directional derivs.
nsuccess = 0;				% nsuccess counts number of successes.
lambda = 1.0;				% Initial scale parameter.
lambdamin = 1.0e-15; 
lambdamax = 1.0e100;
j = 1;					% j counts number of iterations.
xval = x;
vfail = 0;                              % number if failed validations
vx=x;

% Main optimization loop.
while (j <= niters)

  % Calculate first and second directional derivatives.
  if (success == 1)
    mu = d*gradnew';
    if (mu >= 0)
      d = - gradnew;
      mu = d*gradnew';
    end
    kappa = d*d';
    if kappa < eps
      return
    end
    sigma = sigma0/sqrt(kappa);
    xplus = x + sigma*d;
    gplus = feval(gradf, xplus, varargin1{:});
    gamma = (d*(gplus' - gradnew'))/sigma;
  end

  % Increase effective curvature and evaluate step size alpha.
  delta = gamma + lambda*kappa;
  if (delta <= 0) 
    delta = lambda*kappa;
    lambda = lambda - gamma/kappa;
  end
  alpha = - mu/delta;
  
  % Calculate the comparison ratio.
  xnew = x + alpha*d;
  fnew = feval(f, xnew, varargin1{:});
  iter = iter + 1;
  Delta = 2*(fnew - fold)/(alpha*mu);
  if (Delta  >= 0)
    success = 1;
    x = xnew;
    fnow = fnew;
  else
    success = 0;
    fnow = fold;
  end

  % Compute validation value
  vnow=feval(f, x, varargin2{:});
  if vnow<vold
    vold=vnow;
    vx=x;
    vfail=0;
  else
    vfail=vfail+1;
  end
  if display > 1
    fprintf(1, 'Cycle %4d  Error  %8.4f VError %8.4f \n', ...
            j, fold, vnow);
  end
  fs(j)=fnow;
  vs(j)=vnow;
%  subplot(2,1,1)
%  plot(1:j,fs)
%  subplot(2,1,2)
%  plot(1:j,vs)
%  drawnow

  if (success == 1)
    
    % Test for early stop termination
    if vfail>opt.maxfail
      x=vx;
      if (display > 0)
        disp('Early stopping')
      end
      return;
      
      % Test for normal termination
    elseif (max(abs(alpha*d)) < opt.tolx & max(abs(fnew-fold)) < opt.tolfun)
      if (display > 0)
        disp('Tolx and tolfun reached')
      end
      return;

    else
      % Update variables for new position
      fold = fnew;
      gradold = gradnew;
      gradnew = feval(gradf, x, varargin1{:});
      % If the gradient is zero then we are done.
      if (gradnew*gradnew' == 0)
        if (display > 0)
          disp('Gradient zero');
        end
	return;
      end
    end
  end

  % Adjust lambda according to comparison ratio.
  if (Delta < 0.25)
    lambda = min(4.0*lambda, lambdamax);
  end
  if (Delta > 0.75)
    lambda = max(0.5*lambda, lambdamin);
  end

  % Update search direction using Polak-Ribiere formula, or re-start 
  % in direction of negative gradient after nparams steps.
  if (nsuccess == nparams)
    d = -gradnew;
    nsuccess = 0;
  else
    if (success == 1)
%      beta = (gradold - gradnew)*gradnew'/(gradold*gradold');
            beta = (gradold - gradnew)*gradnew'/(mu);
      d = beta*d - gradnew;
    end
  end
  j = j + 1;
end

% If we get here, then we haven't terminated in the given number of 
% iterations.
if (display > 0)
  disp('Warning: Maximum number of iterations has been exceeded');
end
