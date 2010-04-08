function idxs = subv2ind(sz,data)
% SUBV2IND Like the built-in sub2ind, but the subscripts are given as row vectors.
%
% idxs = subv2ind(sz,data)
%
% Copyright (C) 2007, Marcel van Gerven
%

idxs = sum( repmat(cumprod([1 sz(1:(end-1))]),size(data,1),1) .* (data - 1), 2) + 1;

