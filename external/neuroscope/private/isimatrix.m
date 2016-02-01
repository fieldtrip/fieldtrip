%isimatrix - Test if parameter is a matrix of integers (>= 2 columns).
%
%  USAGE
%
%    test = isimatrix(x,test1,test2,...)
%
%    x              parameter to test
%    test1...       optional list of additional tests
%
%  EXAMPLES
%
%    % Test if x is a matrix of doubles
%    isimatrix(x)
%
%    % Test if x is a matrix of strictly positive doubles
%    isimatrix(x,'>0')
%
%  NOTE
%
%    The tests ignore NaNs, e.g. isimatrix([500 nan;4 79]), isimatrix([1 nan 3],'>0') and
%    isimatrix([nan -7;nan nan;-2 -5],'<=0') all return 1.
%
%  SEE ALSO
%
%    See also isdmatrix, isdvector, isdscalar, isivector, isiscalar, isstring,
%    islscalar, islvector, islmatrix.
%

% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = isimatrix(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isimatrix">isimatrix</a>'' for details).');
end

% Test: doubles, two dimensions, two or more columns?
test = isa(x,'double') & length(size(x)) == 2 & size(x,2) >= 2;

% Ignore NaNs (this reshapes the matrix, but it does not matter for the remaining tests)
x = x(~isnan(x));

% Test: integers?
test = test & all(round(x)==x);

% Optional tests
for i = 1:length(varargin),
	try
		if ~eval(['all(x(:)' varargin{i} ');']), test = false; return; end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isimatrix">isimatrix</a>'' for details).']);
	end
end
