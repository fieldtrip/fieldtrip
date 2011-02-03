function [yp,conf] = rls_predict(model, x)
% function yp = rls_predict(model, x)
%
% Predict response from given RLS model and input vector x

% (C) 2010 S.Klanke

yp = model.beta'*x;

if nargout>1
   conf = sqrt(model.mse * (1+x'*model.invX*x));
end

% possible extension:
% also provide predictive variance
% yp +/- sigma*sqrt(1 + x'*model.invX*x)
% sigma would be noise estimate, could be estimated from training data
% and accumulated along with invH and beta. But this only makes sense
% with a forgetting factor (early residuals will be way off).
