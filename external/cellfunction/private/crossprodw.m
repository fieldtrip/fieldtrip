function [y] = crossprodw(x, w)

% expand, compute cross-product and whiten
y = w*crossprod(x);
