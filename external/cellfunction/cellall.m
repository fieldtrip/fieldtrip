function [output] = cellall(x, val)

if nargin==1, val = 1; end
for k = 1:prod(size(x))
  tmp(k) = all(x{k}==val);
end
output = all(tmp==1); 
