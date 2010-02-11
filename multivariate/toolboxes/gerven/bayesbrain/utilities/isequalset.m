function m = isequalset(A,B)
% ISEQUALSET checks if two sets are equal. Note that this only
% works for small positive integers.
%
% m = isequalset(A,B)
%
% Copyright (C) 2007, Marcel van Gerven
%
% $Log: isequalset.m,v $
% Revision 1.1.1.1  2008/02/27 14:42:54  roboos
% Marcel van Gerven, version 27 Feb 2008
%


if ~isequal(max(A),max(B))
    
    m = false;
    
else

    Abits(A) = true;
    Bbits(B) = true;

    m = ~any(xor(Abits,Bbits));    
end