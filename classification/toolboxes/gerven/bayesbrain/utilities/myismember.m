function [x,i] = myismember(A,B)
% ISMEMBER replaces builtin ISMEMBER function. Note that this only
% works for small positive integers and i is the lowest instead of highest
% index!
%
% [x,i] = myismember(A,B)
%
% Copyright (C) 2007, Marcel van Gerven
%
% $Log: myismember.m,v $
% Revision 1.1.1.1  2008/02/27 14:42:54  roboos
% Marcel van Gerven, version 27 Feb 2008
%

if numel(B) == 1

    x = (A == B);
    i = double(x);
else
    x = false(size(A));
    i = zeros(size(A));
    for n=1:numel(A)
        found = find(A(n)==B(:),1);
        if ~isempty(found)
            x(n) = true;
            i(n) = found;
        end
    end
end
