% fastif() - fast if function.
%
% Usage:
%  >> res = fastif(test, s1, s2);
%
% Input:
%   test   - logical test with result 0 or 1
%   s1     - result if 1
%   s2     - result if 0
%
% Output:
%   res    - s1 or s2 depending on the value of the test
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function res = fastif(s1, s2, s3);

if s1
	res = s2;
else
	res = s3;
end;
return;
