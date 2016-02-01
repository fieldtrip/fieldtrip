%isdvector - Test if parameter is a vector of doubles satisfying an optional list of tests.
%
%  USAGE
%
%    test = isdvector(x,test1,test2,...)
%
%    x              parameter to test
%    test1...       optional list of additional tests (see examples below)
%
%  EXAMPLES
%
%    % Test if x is a vector of doubles
%    isdvector(x)
%
%    % Test if x is a vector of strictly positive doubles
%    isdvector(x,'>0')
%
%    % Test if x is a vector of doubles included in [2,3]
%    isdvector(x,'>=2','<=3')
%
%    % Special test: test if x is a vector of doubles of length 3
%    isdvector(x,'#3')
%
%    % Special test: test if x is a vector of strictly ordered doubles
%    isdvector(x,'>')
%
%  NOTE
%
%    The tests ignore NaNs, e.g. isdvector([5e-3 nan]), isdvector([1.7 nan 3],'>0') and
%    isdvector([nan -7.4],'<=0') all return 1.
%
%  SEE ALSO
%
%    See also isdmatrix, isdscalar, isimatrix, isivector, isiscalar, isstring,
%    islscalar, islvector, islmatrix.
%

% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = isdvector(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isdvector">isdvector</a>'' for details).');
end

% Test: double, vector
test = isa(x,'double') & isvector(x);

% Ignore NaNs
x = x(~isnan(x));

% Optional tests
for i = 1:length(varargin),
	try
		if varargin{i}(1) == '#',
			if length(x) ~= str2num(varargin{i}(2:end)), test = false; return; end
		elseif isstring(varargin{i},'>','>=','<','<='),
			dx = diff(x);
			if ~eval(['all(0' varargin{i} 'dx);']), test = false; return; end
		else
			if ~eval(['all(x' varargin{i} ');']), test = false; return; end
		end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isdvector">isdvector</a>'' for details).']);
	end
end
