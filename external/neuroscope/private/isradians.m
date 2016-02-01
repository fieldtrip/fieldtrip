%isradians - Test if parameter is in range [0,2pi] or [-pi,pi].
%
%  USAGE
%
%    test = isradians(x)
%
%    x              array to test (NaNs are ignored)
%
%  OUTPUT
%
%    range          0 if uncertain (issues a warning)
%                   1 for [-pi,pi]
%                   2 for [0,2pi]
%  SEE ALSO
%
%    See also wrap.
%

% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = isradians(x)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help isradians">isradians</a>'' for details).');
end

if ~isa(x,'double'), test = 0; return; end

% Ignore NaN
x = x(~isnan(x));
if isempty(x), test = 0; return; end

% Min and max
x = x(:);
m = min(x);
M = max(x);

% Range test
if m >= -pi && M <= pi,
	test = 1;
elseif m >=0 && M <= 2*pi,
	test = 2;
else
	warning('Angles are neither in [0,2pi] nor in [-pi,pi] (make sure they are in radians).');
	test = 0;
end
