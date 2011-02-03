function crc = neuralynx_crc(dat, dim)

% NEURALYNX_CRC computes a cyclic redundancy check
%
% Use as
%   crc = neuralynx_crc(dat)
%
% Note that the CRC is computed along the first dimension.

% Copyright (C) 2007, Robert Oostenveld
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
% $Id$


nchans   = size(dat,1);
nsamples = size(dat,2);

crc = zeros(1,nsamples,class(dat));

for i=1:nchans
  crc = bitxor(crc, dat(i,:));
end

