function [fr, ut, u1, s] = svdfft(f, n, trltapcnt)

% SVDFFT computes a rotated FFT matrix, using the real part of the cross-spectral
% density matrix. This rotation ensures that the phase relationship of the underlying
% sources does not change, while rotating the channels such that the first channel
% contains the maximal amplitude signal.
%
% Use as
%   [fr, ut] = svdfft(f, n, trltapcnt);
% where
%   n           number of components (orientations) to keep in the output (e.g. 1, 2 or 3)
%   trltapcnt   vector of length Ntrials with the number of tapers
%
% See also SVD

% Copyright (C) 2005-2020, Robert Oostenveld & Jan-Mathijs Schoffelen
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

if nargin == 1
  n         = size(f,1);
  trltapcnt = ones(size(f,1),1);
elseif nargin == 2
  trltapcnt = ones(size(f,1),1);
elseif nargin == 3
  if isempty(n), n = size(f,1); end
end

if all(trltapcnt==trltapcnt(1))
  c = f * f';
else
  trltapcnt = trltapcnt(:);
  sumtapcnt = cumsum([0;trltapcnt]);
  c         = zeros(size(f,1), size(f,1));
  for j = 1:length(sumtapcnt)-1
    ftap = f(:, sumtapcnt(j)+1:sumtapcnt(j+1));
    c = c + ftap * ftap' ./ trltapcnt(j);
  end
end

if n==size(f,1)
  % do a complete decomposition
  [u, s] = svd(real(c));
elseif n<1
  % do a complete decomposition and only keep the biggest components which together explain n percent of the variance
  [u, s] = svd(real(c));
  s      = cumsum(diag(s))./sum(diag(s));
  n      = length(find(s<=n));
  u      = u(:, 1:n);
else
  % only decompose the first n components
  [u, s] = svds(real(c),n);
end

ut  = u';      % this rotates the data in the direction of the maximum power
fr  = ut * f;  % apply the rotation on the data
u1  = u(:,1);
s   = diag(s);
