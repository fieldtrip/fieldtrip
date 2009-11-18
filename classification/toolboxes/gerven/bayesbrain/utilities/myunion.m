function m = myunion(A,B)
% MYUNION replaces builtin UNION function. Note that this only
% works for small positive integers.
%
% m = myunion(A,B)
%
% Copyright (C) 2007, Marcel van Gerven
%
% $Log: myunion.m,v $
% Revision 1.1.1.1  2008/02/27 14:42:54  roboos
% Marcel van Gerven, version 27 Feb 2008
%

if isempty(A) 
    
    m = B;
    
elseif isempty(B)
    
    m = A;
    
else

    bits = false(1,max(max(A),max(B)));

    bits(A) = true;
    bits(B) = true;
    
    m = find(bits);
end