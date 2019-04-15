function [x, fs, ps] = scg2(f, x, opt, gradf, varargin)
%SCG2	Scaled conjugate gradient optimization
%
%  Description
%    X = SCG2(F, X, OPT, GRADF, P1,P2,...) uses a scaled conjugate
%    gradients algorithm to find a local minimum of the function
%    F(X,P1,P2,...) whose gradient is given by GRADF(X,P1,P2,...).
%
%    Here X is a row vector and F returns a scalar value. The point
%    at which F has a local minimum is returned as X.
%
%    See SCG2_OPT for the optional parameters in the OPT structure.
%
%  See also SCG2_OPT, SCGES

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

%  Check gradients
if opt.checkgrad
  feval('gradcheck', x, f, gradf, varargin{:});
end

sigma0 = 1.0e-4;
fold = feval(f, x, varargin{:});	% Initial function value.
iter = 0;
gradnew = feval(gradf, x, varargin{:}); % Initial gradient.
gradold = gradnew;
d = - gradnew;				% Initial search direction.
success = 1;                            % Force calculation of directional derivs.
nsuccess = 0;				% nsuccess counts number of successes.
lambda = 1.0;				% Initial scale parameter.
lambdamin = 1.0e-15; 
lambdamax = 1.0e100;
j = 1;					% j counts number of iterations.
if nargout >= 2
  fs(j, :) = fold;
  if nargout == 3
    ps(j, :) = x;
  end
end

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
    gplus = feval(gradf, xplus, varargin{:});
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
  fnew = feval(f, xnew, varargin{:});
  iter = iter + 1;
  Delta = 2*(fnew - fold)/(alpha*mu);
  if (Delta  >= 0)
    success = 1;
    nsuccess = nsuccess + 1;
    x = xnew;
    fnow = fnew;
  else
    success = 0;
    fnow = fold;
  end

  if nargout >= 2
    % Store relevant variables
    fs(j) = fnow;		% Current function value
    if nargout >= 3
      ps(j,:) = x;	% Current position
    end
  end    
  if display > 1
    fprintf(1, 'Cycle %4d  Error %8.4f  Scale %8.4f\n', j, fnow, lambda);
  end
  
  if (success == 1)
    
    % Test for termination
    if (max(abs(alpha*d)) < opt.tolx & max(abs(fnew-fold)) < opt.tolfun)
      if (display > 0)
        disp('Tolx and tolfun reached')
      end
      return;

    else
      % Update variables for new position
      fold = fnew;
      gradold = gradnew;
      gradnew = feval(gradf, x, varargin{:});
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
