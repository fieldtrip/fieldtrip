function x = trimnumericalerror(x, tolerance, replace)

% TRIMNUMERICALERROR replaces values in the input array that are likely to result in
% numerical errors. It looks at the ratio of all values to the largest value, values
% smaller than the relative tolerance are replaced.

if nargin<2
  tolerance = eps * 1e3;
end

if nargin<3
  replace = 0;
end

sel = abs(x) < tolerance * max(abs(x));
x(sel) = replace;
