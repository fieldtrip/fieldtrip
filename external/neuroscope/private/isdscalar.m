%isdscalar - Test if parameter is a scalar (double) satisfying an optional list of tests.
%
%  USAGE
%
%    test = isdscalar(x,test1,test2,...)
%
%    x              parameter to test
%    test1...       optional list of additional tests
%
%  EXAMPLES
%
%    % Test if x is a scalar (double)
%    isdscalar(x)
%
%    % Test if x is a strictly positive scalar (double)
%    isdscalar(x,'>0')
%
%    % Test if x is a scalar (double) included in [2,3]
%    isdscalar(x,'>=2','<=3')
%
%  NOTE
%
%    The tests ignore NaN, e.g. isdscalar(nan), isdscalar(nan,'>0') and isdscalar(nan,'<=0')
%    all return 1.
%
%  SEE ALSO
%
%    See also isdmatrix, isdvector, isimatrix, isivector, isiscalar, isstring,
%    islscalar, islvector, islmatrix.
%

% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = isdscalar(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isdscalar">isdscalar</a>'' for details).');
end

% Test: double, scalar
test = isa(x,'double') & isscalar(x);
if ~test, return; end

% Optional tests
for i = 1:length(varargin),
	try
		if ~eval(['x' varargin{i} ';']), test = false; return; end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isdscalar">isdscalar</a>'' for details).']);
	end
end
