% NANSUM provides a replacement for MATLAB's nanmean.
%
% For usage see SUM.

function y = nansum(x, dim)

if nargin == 1
  dim = 1;
end
idx = isnan(x);
x(idx) = 0;
y = sum(x, dim);
end