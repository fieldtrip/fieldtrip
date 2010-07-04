function d = idummy(x, method)
%DUMMY  Creates a matrix of dummy variables for a discrete variable
%   D=IDUMMY(X,METHOD) creates an array D of dummy variables for the
%   grouping variable I (integers 1,...,g), using the method specified:
%
%   method = 1:   0/-1/1 coding, full rank
%   method = 2:   0/1 coding, full rank
%   method = 3:   0/1 coding, overdetermined

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.40 $  $Date: 2009/01/23 15:13:49 $

if (nargin < 2)
   method = 1;
end

n = length(x);
g = max(x);
ncols = g - (method ~= 3);
d = repmat(0, n, ncols);

if (g > 1)
   % Fill in -1 for the first level
   if (method == 1)
      i = find(x == 1);
      d(find(x == 1),:) = -1;
   end
   
   % Fill in 1 in the appropriate column for other levels
   m3 = (method == 3);
   for j=(2-m3):g
      d(find(x == j),j-1+m3) = 1;
   end
end
