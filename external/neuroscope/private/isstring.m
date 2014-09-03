%isstring - Test if parameter is an (admissible) character string.
%
%  USAGE
%
%    test = isstring(x,string1,string2,...)
%
%    x              item to test
%    string1...     optional list of admissible strings
%
%  SEE ALSO
%
%    See also isdmatrix, isdvector, isdscalar, isimatrix, isivector, isiscalar.
%

% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = isstring(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isstring">isstring</a>'' for details).');
end

test = true;

if ~ischar(x),
	test = false;
	return;
end

if isempty(varargin), return; end

for i = 1:length(varargin),
	if strcmp(x,varargin{i}), return; end
end

test = false;
