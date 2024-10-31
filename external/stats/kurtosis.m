function k = kurtosis(x, flag, dim)

% K = KURTOSIS(X, FLAG, DIM) computes the kurtosis of the values in X along the 
% dimention DIM. FLAG determines whether the bias is corrected. (0: bias
% is corrected, 1: bias is not corrected (default)).

if nargin < 2 || isempty(flag)
  flag = 1;
end

if nargin < 3 || isempty(dim)
  dim = find(size(x)>1,1,'first');
end

% Center X, compute its fourth and second moments, and compute the
% uncorrected kurtosis.
x0 = x - mean(x, dim);
s2 = mean(x0.^2, dim); % this is the biased variance estimator
m4 = mean(x0.^4, dim);
k = m4 ./ s2.^2;

% Bias correct the kurtosis.
if ~flag
  n = sum(~isfinite(x), dim);
  n(n<4) = NaN; % bias correction is not defined for n < 4.
  k = ((n+1).*k - 3.*(n-1)) .* (n-1)./((n-2).*(n-3)) + 3;
end
