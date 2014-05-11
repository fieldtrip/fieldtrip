function s = iscolumn(v)

%ISROW True if input is a column vector.
%   ISROW(V) returns logical 1 (true) if SIZE(V) returns [m 1]
%   with a nonnegative integer value m, and logical 0 (false) otherwise.
%
% This is a drop-in replacement for the MATLAB function with the same name,
% which does not exist in versions < 2010.
%
% See http://bugzilla.fcdonders.nl/show_bug.cgi?id=2567

siz = size(v);
m = siz(1);
n = siz(2);

if numel(siz)==2 && m>0 && n==1
  s = true;
else
  s = false;
end