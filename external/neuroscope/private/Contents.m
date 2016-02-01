% Helper functions for FMAToolbox.
%
% Parameter testing.
%
%   isstring              - Test if parameter is an (admissible) character string.
%   isdscalar             - Test if parameter is a scalar (double) satisfying an optional list of tests.
%   isdvector             - Test if parameter is a vector of doubles satisfying an optional list of tests.
%   isdmatrix             - Test if parameter is a matrix of doubles (>= 2 columns).
%   isiscalar             - Test if parameter is a scalar (integer) satisfying an optional list of tests.
%   isivector             - Test if parameter is a vector of integers satisfying an optional list of tests.
%   isimatrix             - Test if parameter is a matrix of integers (>= 2 columns).
%   islscalar             - Test if parameter is a (pseudo) logical scalar.
%   islvector             - Test if parameter is a (pseudo) logical vector satisfying an optional list of tests.
%   islmatrix             - Test if parameter is a (pseudo) logical matrix (>= 2 columns).
%   isradians             - Test if parameter is in range [0,2pi] or [-pi,pi].
%   wrap                  - Set radian angles in range [0,2pi] or [-pi,pi].
%