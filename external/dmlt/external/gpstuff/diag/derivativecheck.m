function delta  = derivativecheck(w, fun, epsilon)
%DERIVATIVECHECK Compare user-supplied derivatives to
%                finite-differencing derivatives.
%
%  Description
%    This function is intended as a utility to check whether a
%    gradient calculation has been correctly implemented for a
%    given function. 
%
%    DERIVATIVECHECK(X, FUN) checks how accurate the user-supplied
%    derivatives of the function FUN are at X. FUN accepts input X
%    and returns a scalar function value F and its scalar or vector
%    gradient G evaluated at X A central difference formula with
%    step size 1.0e-6 is used, and the results for both gradient
%    function and finite difference approximation are printed.
%
%    DERIVATIVECHECK(X, FUN, EPSILON) computes the finite
%    difference gradients using user given step size EPSILON.
%
%    DELTA=DERIVATIVECHECK(X, FUN) returns the delta between the
%    user-supllied and finite difference derivatives.
%

% Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)
% Copyright (c) 2010-2011 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% Reasonable value for step size
if nargin<3
  epsilon = 1.0e-5;
end

% Treat
nparams = length(w);
%deltaf = zeros(1, nparams);
step0 = zeros(1, nparams);

for i = 1:nparams
  % Move a small way in the ith coordinate of w
  step = step0;
  step(i) = 1;
  fplus  = fun(w+epsilon.*step);
  fminus = fun(w-epsilon.*step);
  % Use central difference formula for approximation
  if nparams>1 && numel(fplus)==1
    deltaf(1,i) = 0.5*(fplus - fminus)/epsilon;
  elseif nparams>1 && numel(fplus)>1
    deltaf(:,i) = 0.5*(fplus(i)' - fminus(i)')/epsilon;
  else
    deltaf(i,:) = 0.5*(fplus(:)' - fminus(:)')/epsilon;
  end
end
[tmp,gradient] = fun(w);
gradient=gradient(:)';
fprintf(1, 'Checking gradient ...\n\n');
fprintf(1, '   analytic   diffs     delta\n\n');
disp([gradient', deltaf', gradient' - deltaf'])
delta = gradient' - deltaf';
