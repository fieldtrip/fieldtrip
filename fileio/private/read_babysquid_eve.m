function [smp, tim, val3, val4] = read_babysquid_eve(filename)

% READ_BABYSQUID_EVE imports events from the *.eve file that accompanies the *.fif
% file.
%
% Use as
%  [smp, tim, val3, val4] = read_babysquid_eve(filename)

% Copyright (C) 2013 Robert Oostenveld
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

% TODO check sample indexing (1 or 0 based)
% TODO clarify val3 and val4

fid = fopen(filename, 'rt');
[a, count] = fscanf(fid, '%f', [inf]);
fclose(fid);

% the output of fscanf can be reshaped into a transposed matrix
m = 4;
n = count/m;
a = reshape(a, m, n);

smp  = a(1,:)'; % sample indices should start counting at 1, not at 0
tim  = a(2,:)';
val3 = a(3,:)';
val4 = a(4,:)';

