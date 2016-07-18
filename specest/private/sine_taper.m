function d = sine_taper(n, k)

% Compute Riedel & Sidorenko sine tapers.
% sine_taper(n, k) produces the first 2*k tapers of length n,
% returned as the columns of d.

% Copyright (C) 2006, Tom Holroyd
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin < 2
  error('usage: sine_taper(n, k)');
end

k = round(k * 2);
if k <= 0 || k > n
  error('sine_taper: k is %g, must be in (1:n)/2', k)
end

x = (1:k) .* (pi / (n + 1));
d = sqrt(2 / (n + 1)) .* sin((1:n)' * x);

return
