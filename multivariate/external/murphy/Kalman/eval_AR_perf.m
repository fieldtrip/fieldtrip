function [ypred, ll, mse] = eval_AR_perf(coef, C, y, model)
% Evaluate the performance of an AR model.
% 
% Inputs
% coef(:,:,k,m) - coef. matrix to use for k steps back, model m
% C(:,:,m)      - cov. matrix for model m
% y(:,t)        - observation at time t
% model(t)      - which model to use at time t (defaults to 1 if not specified)
%
% Outputs
% ypred(:,t)    - the predicted value of y at t based on the evidence thru t-1.
% ll            - log likelihood
% mse           - mean squared error = sum_t d_t . d_t, where d_t = pred(y_t) - y(t)

[s T] = size(y);
k = size(coef, 3);
M = size(coef, 4);

if nargin<4, model = ones(1, T); end

ypred = zeros(s, T);
ypred(:, 1:k) = y(:, 1:k);
mse = 0;
ll = 0;
for j=1:M
  c(j) = log(normal_coef(C(:,:,j)));
  invC(:,:,j) = inv(C(:,:,j));
end
coef = reshape(coef, [s s*k M]);

for t=k+1:T
  m = model(t-k);
  past = y(:,t-1:-1:t-k);
  ypred(:,t) = coef(:, :, m) * past(:);
  d = ypred(:,t) - y(:,t);
  mse = mse + d' * d;
  ll = ll + c(m) - 0.5*(d' * invC(:,:,m) * d);
end
mse = mse / (T-k+1);

