function m = mysetdiff(A,B)
% MYSETDIFF replaces builtin SETDIFF function. Note that this only
% works for small positive integers.
%
% m = mysetdiff(A,B)
%
% Copyright (C) 2007, Marcel van Gerven
%
% $Log: mysetdiff.m,v $
% Revision 1.1.1.1  2008/02/27 14:42:54  roboos
% Marcel van Gerven, version 27 Feb 2008
%

if isempty(B)
    
    m = A;
    
elseif isempty(A)
    
    m = [];

else

    bits = true(1,max(max(A),max(B)));

    bits(B) = false;
    m = A(bits(A));
end