function model = rls_update(model, x, y)
% function model = rls_update(model, x, y)
%
% Update a recursive least squares model with a new training tuple (x,y)
% which must both be column vectors.

% (C) 2010 S.Klanke

% Simple but slow equations:
% model.H    = model.H + x*x';
% model.XtY  = model.XtY  + x*y';
% model.beta = inv(model.H) * model.XtY;
    
invHx = model.invH*x;

L = invHx/(model.lambda + x'*invHx);

model.invH = (model.invH - L*invHx')/model.lambda;
model.beta = model.beta + L*( y' - x'*model.beta );

N = 0.995*model.N; % slightly forget old errors

err = (model.beta'*x - y).^2;
model.mse = (N*model.mse + err)/(N+1);
model.N = N+1;

% possible extension:
% use a forgetting factor lambda instead of "1" in the computation of L
% lambda < 1 downweights old samples (exponential forgetting)



% RLS is based on the Woodbury identity (see Matrix Cookbook)
% inv(A + BC) = inv(A) - inv(A)*B*inv(eye + C*inv(A)*B)*C*inv(A)
%  A = H, B = x,  C = x'
% inv(H + xx') = invH - invH*x*(1 + x'*invH*x)^(-1)*x'*invH'
%              = invH - invHx*(1+x'invHx)^(-1)*invHx'
% 
% and further
%
% beta = inv(H + xx')*(accXtY + x*y') = invH*accXtY + invH*x*y'
% -invHx*()*x'*invH'*accXtY - invHx*()*invHx'*x*y'
% = old_beta + invHx*y' - invHx*()*x'*old_beta 
%   - (x'invHx)/(1+x'invHx)*invHx*y'
% = old_beta + invHx/(1+x'invHx)*(y' - x'*old_beta)
