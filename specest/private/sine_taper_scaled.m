function d = sine_taper_scaled(n, k)

% Compute Riedel & Sidorenko sine tapers.
% sine_taper_scaled(n, k) produces the first 2*k tapers of length n,
% returned as the columns of d. The norm of the tapers will not be 1. The
% norm is a function of the number of the taper in the sequence. This is to
% mimick behavior of the scaling of the resulting powerspectra prior to
% april 29, 2011. Before april 29, 2011, equivalent scaling was applied to
% the powerspectra of the tapered data segments, prior to averaging.

% Copyright (C) 2006, Tom Holroyd
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
% $Id: sine_taper_scaled.m 2885 2011-02-16 09:41:58Z roboos $

if nargin < 2
  error('usage: sine_taper_scaled(n, k)');
end

k = round(k * 2);
if k <= 0 || k > n
  error('sine_taper_scaled: k is %g, must be in (1:n)/2', k)
end

x = (1:k) .* (pi / (n + 1));
d = sqrt(2 / (n + 1)) .* sin((1:n)' * x);

scale = zeros(k);
for i = 1:k
  scale(i,i) = sqrt(1 - ((i-1)./(k-1)).^2);
end
d = d * scale;

return
