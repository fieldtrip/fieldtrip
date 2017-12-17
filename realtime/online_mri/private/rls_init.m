function model = rls_init(dimX, dimY, gamma, lambda)
% function model = rls_init(dimX, dimY [, gamma = 1e-6 [, lambda = 1]])
%
% Initialise a recursive least squares model with input dimension dimX and output dimension dimY.
% The optional third parameter sets the regularisation coefficient (as in ridge regression).
% The optional fourth argument sets the forgetting factor (1 = no forgetting).
%
% See also RLS_UPDATE, RLS_PREDICT

% Straightforward but slow equations:
%
% model.H    = model.H + x*x';
% model.XtY  = model.XtY  + x*y';
% model.beta = inv(model.H) * model.XtY;
%
% RLS works by directly accumulating the inverse of H and also updating beta without 
% computing XtY
%
% (C) 2010 S.Klanke

if nargin<3
  gamma = 1e-6;
elseif gamma <= 0
  error 'Regularisation parameter must be positive'
end

if nargin<4
  lambda = 1;
elseif lambda < 0 || lambda > 1
  error 'Forgetting factor must be within [0;1]';
end
    
model = struct('beta', zeros(dimX,dimY), 'invH', (1/gamma)*eye(dimX), 'N', 0, 'mse', zeros(dimY,1), 'lambda', lambda);
