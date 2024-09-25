function [pos, tri] = read_vtk(fn)

% READ_VTK reads a triangulation from a VTK (Visualisation ToolKit) format file
% Supported are triangles and other polygons.
%
% Use as
%   [pnt, tri] = read_vtk(filename)
%
% See also WRITE_VTK

% Copyright (C) 2002-2024, Robert Oostenveld
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

fid = fopen_or_error(fn, 'rt');

npos = 0;
while (~npos)
  line = fgetl(fid);
  if contains(line, 'POINTS')
    npos = sscanf(line, 'POINTS %d float');
  end
end
pos = zeros(npos, 3);
for i=1:npos
  pos(i,:) = fscanf(fid, '%f', 3)';
end

ntri = 0;
while (~ntri)
  line = fgetl(fid);
  if contains(line, 'POLYGONS')
    tmp = sscanf(line, 'POLYGONS %d %d');
    ntri  = tmp(1);          % number of triangles
    nvert = tmp(2)/ntri - 1; % number of vertices per polygon
  end
end

tri = zeros(ntri, nvert+1);
for i=1:ntri
  tri(i,:) = fscanf(fid, '%d', nvert+1)';
end
% drop the first column
tri = tri(:,2:(nvert+1)) + 1;

fclose(fid);
