function [elec] = read_polhemus_pos(filename)

% READ_POLHEMUS_POS reads electrode positions measured with the Polhemus tracker in
% one of the EEG labs at the DCCN. The software used with the Polhemus is from CTF.
%
% Use as:
%   [elec] = read_polhemus_pos(filename)
%
% This returns an electrode structure with
%   elec.label     cell-array with electrode labels (strings)
%   elec.pnt       position of each electrode

% Copyright (C) 2004-2019, Robert Oostenveld
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

fid = fopen_or_error(filename, 'rt');
line = fgetl(fid);
Nchan = str2double(line);

for i=1:Nchan
  line = fgetl(fid);
  [t, r] = strtok(line);
  elec.label{i} = char(t);
  elec.pnt(i,:) = sscanf(r, '%f')';
end
elec.label = elec.label(:);

try
  % read the fiducials
  line = fgetl(fid);
  [t, r] = strtok(line);
  fiducial.label{1} = char(t);
  fiducial.pnt(1,:) = sscanf(r, '%f')';
  line = fgetl(fid);
  [t, r] = strtok(line);
  fiducial.label{2} = char(t);
  fiducial.pnt(2,:) = sscanf(r, '%f')';
  line = fgetl(fid);
  [t, r] = strtok(line);
  fiducial.label{3} = char(t);
  fiducial.pnt(3,:) = sscanf(r, '%f')';
  % add the fiducials to the electrode array
  elec.label = cat(1, elec.label(:), fiducial.label(:));
  elec.pnt   = cat(1, elec.pnt, fiducial.pnt);
catch
  % do nothing
end

fclose(fid);
