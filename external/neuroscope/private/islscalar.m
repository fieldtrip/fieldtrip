%islscalar - Test if parameter is a (pseudo) logical scalar.
%
%  USAGE
%
%    test = islscalar(x)
%
%    x              parameter to test
%
%  NOTE
%
%    To be considered logical, the scalar should be equal to 0 or 1, but it does
%    not need to actually be of class 'logical' (class(x) could be e.g. 'double').
%
%  SEE ALSO
%
%    See also islvector, islmatrix, isdmatrix, isdvector, isdscalar, isimatrix, isivector,
%    isstring.
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = islscalar(x)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help islscalar">islscalar</a>'' for details).');
end

% Test: double, scalar
test = islogical(x) & isscalar(x);


function test = islogical(x)

test = builtin('islogical',x) | all(x(:)==0|x(:)==1);
