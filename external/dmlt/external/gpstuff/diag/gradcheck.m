function delta  = gradcheck(w, func, grad, varargin)
%GRADCHECK Checks a user-defined gradient function using finite differences.
%
%	Description
%	This function is intended as a utility to check whether a gradient
%       calculation has been correctly implemented for a given function. 
%       GRADCHECK(W, FUNC, GRAD) checks how accurate the gradient GRAD of 
%       a function FUNC is at a parameter vector X.   A central difference 
%       formula with step size 1.0e-6 is used, and the results for both 
%       gradient function and finite difference approximation are printed.
%
%	GRADCHECK(X, FUNC, GRAD, P1, P2, ...) allows additional arguments to
%	be passed to FUNC and GRAD.
%

%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Reasonable value for step size
epsilon = 1.0e-6;

func = fcnchk(func, length(varargin));
grad = fcnchk(grad, length(varargin));

% Treat
nparams = length(w);
deltaf = zeros(1, nparams);
step = zeros(1, nparams);
for i = 1:nparams
  % Move a small way in the ith coordinate of w
  step(i) = 1.0;
  fplus  = feval('linef', epsilon, func, w, step, varargin{:});
  fminus = feval('linef', -epsilon, func, w, step, varargin{:});
  % Use central difference formula for approximation
  deltaf(i) = 0.5*(fplus - fminus)/epsilon;
  step(i) = 0.0;
end
gradient = feval(grad, w, varargin{:});
fprintf(1, 'Checking gradient ...\n\n');
fprintf(1, '   analytic   diffs     delta\n\n');
disp([gradient', deltaf', gradient' - deltaf'])

delta = gradient' - deltaf';

function y = linef(lambda, fn, x, d, varargin)
%LINEF	Calculate function value along a line.
%
%	Description
%	LINEF(LAMBDA, FN, X, D) calculates the value of the function FN at
%	the point X+LAMBDA*D.  Here X is a row vector and LAMBDA is a scalar.
%
%	LINEF(LAMBDA, FN, X, D, P1, P2, ...) allows additional arguments to
%	be passed to FN().   This function is used for convenience in some of
%	the optimization routines.
%
%	See also
%	GRADCHEK, LINEMIN
%

%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Check function string
fn = fcnchk(fn, length(varargin));

y = feval(fn, x+lambda.*d, varargin{:});
