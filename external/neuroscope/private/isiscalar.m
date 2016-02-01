%isiscalar - Test if parameter is a scalar (integer) satisfying an optional list of tests.
%
%  USAGE
%
%    test = isiscalar(x,test1,test2,...)
%
%    x              parameter to test
%    test1...       optional list of additional tests
%
%  EXAMPLES
%
%    % Test if x is a scalar (double)
%    isiscalar(x)
%
%    % Test if x is a strictly positive scalar (double)
%    isiscalar(x,'>0')
%
%    % Test if x is a scalar (double) included in [2,3]
%    isiscalar(x,'>=2','<=3')
%
%  NOTE
%
%    The tests ignore NaN, e.g. isiscalar(nan), isiscalar(nan,'>0') and isiscalar(nan,'<=0')
%    all return 1.
%
%  SEE ALSO
%
%    See also isdmatrix, isdvector, isdscalar, isimatrix, isivector, isstring,
%    islscalar, islvector, islmatrix.
%


% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = isiscalar(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isiscalar">isiscalar</a>'' for details).');
end

% Test: double, scalar
test = isa(x,'double') & isscalar(x);
if ~test, return; end

% Test: integers?
test = test & round(x)==x;

% Optional tests
for i = 1:length(varargin),
	try
		if ~eval(['x' varargin{i} ';']), test = false; return; end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isiscalar">isiscalar</a>'' for details).']);
	end
end
