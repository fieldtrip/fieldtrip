function g = logist(x)

% x: 1 x neval
% g: 1 x neval

g = zeros(size(x));

ii = find(x > 20);
g(ii) = -exp(-x(ii));

ii = find(x < -50);
g(ii) = x(ii);

ii = find(x <=20 & x >= -50);
g(ii) = -log(1 + exp(-x(ii)));