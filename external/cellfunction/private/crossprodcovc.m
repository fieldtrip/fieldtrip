function [c,m] = crossprodcovc(x, p)

x     = crossprod(x, p);
[c,m] = covc(x,2);