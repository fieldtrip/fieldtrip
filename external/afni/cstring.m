function str = cstring(x)
%
%CSTRING Print to a string like c
%  CSTRING(x) takes an array x, of type char, and prints it to
%  a string the way C would.
%  the first null, or non-ascii character is the terminator

i=1;
while (x(i)>0 & x(i)<128 & i<length(x))
  i=i+1;
end

if i==1
  str = char(0);
else
  str = char(x(1:i-1));
end

str = sprintf(str,'%s');

return
