function ab = termcross(a,b)
%TERMCROSS Multiply dummy variables for two terms to get interaction

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.40 $  $Date: 2009/01/23 15:13:50 $
if (isempty(a)), ab = b; return, end
if (isempty(b)), ab = a; return, end

na = size(a,2);
nb = size(b,2);
acols = repmat((1:na), 1, nb);
bcols = reshape(repmat((1:nb), na, 1), 1, na*nb);
ab = a(:,acols) .* b(:,bcols);
