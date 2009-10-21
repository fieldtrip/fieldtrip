function d = sine_taper(n, k)

% Compute Riedel & Sidorenko sine tapers.
% sine_taper(n, k) produces the first 2*k tapers of length n,
% returned as the columns of d.

% Copyright (C) 2006, Tom Holroyd
%
% $Log: sine_taper.m,v $
% Revision 1.1  2006/03/06 09:28:38  roboos
% new implementation by Tom Holroyd
%

if nargin < 2
  error('usage: sine_taper(n, k)');
end

k = round(k * 2);
if k <= 0 || k > n
  error('sine_taper: k is %g, must be in (1:n)/2', k)
end

x = (1:k) .* (pi / (n + 1));
d = sqrt(2 / (n + 1)) .* sin((1:n)' * x);

return
