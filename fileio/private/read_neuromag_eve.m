function [col1, col2, col3, col4] = read_neuromag_eve(filename)

% READ_NEUROMAG_EVE imports events from the *.eve marker file that can accompany a
% *.fif dataset.
%
% Use as
%  [smp, tim, val3, val4] = read_neuromag_eve(filename)
%
% Column one is the sample number. Column two is the time. Column three is is most
% cases always zero, but is useful when you need to mark a segment rather than a
% time point. Column four value is the event type you assign, i.e. the value of
% the trigger.
%
% The recording of the data to disk may start later than the actual data
% acquisition. This is represented in hdr.orig.raw.first_samp. This potential
% offset needs to be taken into acocunt when combining it with the data from the
% file on disk.

% Copyright (C) 2013, Robert Oostenveld
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

% TODO check sample indexing (1 or 0 based)

fid = fopen_or_error(filename, 'rt');
[a, count] = fscanf(fid, '%f', [inf]);
fclose(fid);

% the output of fscanf can be reshaped into a transposed matrix
m = 4;
n = count/m;
a = reshape(a, m, n);

col1 = a(1,:)'; % sample indices start counting at 0
col2 = a(2,:)'; % time is not guaranteed to start with 0 at the start of the file
col3 = a(3,:)';
col4 = a(4,:)';
