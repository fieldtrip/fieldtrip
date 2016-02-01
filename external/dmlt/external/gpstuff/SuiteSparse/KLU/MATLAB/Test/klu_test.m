function klu_test (nmat)
%klu_test KLU test
% Example:
%   klu_test
%
% See also klu

% Copyright 2004-2007 Timothy A. Davis, Univ. of Florida
% http://www.cise.ufl.edu/research/sparse

if (nargin < 1)
    nmat = 500 ;
end

test1 (nmat) ;
test2 (nmat) ;
test3 ;
test4 (nmat) ;
test5  ;
