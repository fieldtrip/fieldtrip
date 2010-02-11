function idxs = subv2ind(sz,data)
% SUBV2IND Like the built-in sub2ind, but the subscripts are given as row vectors.
%
% idxs = subv2ind(sz,data)
%
% Copyright (C) 2007, Marcel van Gerven
%
% $Log: subv2ind.m,v $
% Revision 1.1.1.1  2008/02/27 14:42:54  roboos
% Marcel van Gerven, version 27 Feb 2008
%

idxs = sum( repmat(cumprod([1 sz(1:(end-1))]),size(data,1),1) .* (data - 1), 2) + 1;

