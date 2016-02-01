%islvector - Test if parameter is a (pseudo) logical vector satisfying an optional list of tests.
%
%  USAGE
%
%    test = islvector(x,test1,test2,...)
%
%    x              parameter to test
%    test1...       optional list of additional tests (see examples below)
%
%  EXAMPLES
%
%    % Test if x is a logical vector
%    islvector(x)
%
%    % Special test: test if x is a logical vector of length 3
%    islvector(x,'#3')
%
%  NOTE
%
%    To be considered logical, the vector should contain only values 0 and 1, but it
%    does not need to actually be of class 'logical' (class(x) could be e.g. 'double').
%
%  SEE ALSO
%
%    See also islscalar, islmatrix, isdmatrix, isdvector, isdscalar, isimatrix, isiscalar,
%    isstring.
%

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = islvector(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help islvector">islvector</a>'' for details).');
end

% Test: logical, vector
test = islogical(x) & isvector(x);

% Optional tests
for i = 1:length(varargin),
	try
		if varargin{i}(1) == '#',
			if length(x) ~= str2num(varargin{i}(2:end)), test = false; return; end
		end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help islvector">islvector</a>'' for details).']);
	end
end

function test = islogical(x)

test = builtin('islogical',x) | all(x(:)==0|x(:)==1);
