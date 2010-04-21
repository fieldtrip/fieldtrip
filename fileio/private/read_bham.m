function [dat, lab] = read_bham(filename)

% READ_BHAM reads the EEG data files as recorded by Praamstra in Birmingham
% the datafiles are in a particular ascii format
%
% [dat, lab] = read_bham(filename)

% Copyright (C) 2000, Robert Oostenveld
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

fid = fopen(filename, 'rt');

lablen = 6;
line   = fgetl(fid);
numelc = 0;
while ~isempty(line)
  numelc = numelc + 1;
  [t, r] = strtok(line);
  lab(numelc,:) = [blanks(lablen-length(t)), t];
  line = r;
end

buf = fscanf(fid, '%f');
dat = zeros(numelc, length(buf)/numelc);
dat(:) = buf;

