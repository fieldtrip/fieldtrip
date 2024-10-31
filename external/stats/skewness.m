function s = skewness(x, flag, dim)

% S = SKEWNESS(X, FLAG, DIM) computes the skewness of the values in X along the 
% dimention DIM. FLAG determines whether the bias is corrected. (0: bias
% is corrected, 1: bias is not corrected (default)).

if nargin < 2 || isempty(flag)
  flag = 1;
end

if nargin < 3 || isempty(dim)
  dim = find(size(x)>1,1,'first');
end

% Center X, compute its third and second moments, and compute the
% uncorrected skewness.
% x0 = x - repmat(nanmean(x,dim), t);
x0 = x - mean(x, dim, 'omitnan');
s2 = mean(x0.^2, dim, 'omitnan'); % this is the biased variance estimator
m3 = mean(x0.^3, dim, 'omitnan');
s  = m3 ./ s2.^(1.5);
% Bias correct the skewness.
if ~flag
  n = sum(~isfinite(x), dim);
  n(n<3) = NaN; % bias correction is not defined for n < 3.
  s = s .* sqrt((n-1)./n) .* n./(n-2);
end
