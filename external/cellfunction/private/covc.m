function [c,mx,my] = covc(x, y, dim)

if nargin==2,
  dim = y;
  y   = x;
end

if dim==1,
  c  = x'*y;
elseif dim==2,
  c = x*y';
end

mx = sum(x,dim);
my = sum(y,dim);
