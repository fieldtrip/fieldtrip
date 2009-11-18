function m = myintersect(A,B)
% MYINTERSECT replaces builtin INTERSECT function. Note that this only
% works for small positive integers.
%
% m = myintersect(A,B)
%
% Copyright (C) 2007, Marcel van Gerven
%
% $Log: myintersect.m,v $
% Revision 1.1.1.1  2008/02/27 14:42:54  roboos
% Marcel van Gerven, version 27 Feb 2008
%

if isempty(A) || isempty(B)
    
    m = [];
    
else

    bits = false(1,max(max(A),max(B)));

    bits(B) = true;
    m = A(bits(A));
end