function [pnt, tri] = read_vtk(fn)

% READ_VTK reads a triangulation from a VTK (Visualisation ToolKit) format file
% Supported are triangles.
%
% Use as
%   [pnt, tri] = read_vtk(filename)
%
% See also WRITE_VTK

% Copyright (C) 2002, Robert Oostenveld
%
% $Id$

fid = fopen_or_error(fn, 'rt');

npnt = 0;
while (~npnt)
  line = fgetl(fid);
  if contains(line, 'POINTS')
    npnt = sscanf(line, 'POINTS %d float');
  end
end
pnt = zeros(npnt, 3);
for i=1:npnt
  pnt(i,:) = fscanf(fid, '%f', 3)';
end

ntri = 0;
while (~ntri)
  line = fgetl(fid);
  if contains(line, 'POLYGONS')
    tmp = sscanf(line, 'POLYGONS %d %d');
    ntri = tmp(1);
  end
end
tri = zeros(ntri, 4);
for i=1:ntri
  tri(i,:) = fscanf(fid, '%d', 4)';
end
tri = tri(:,2:4) + 1;

fclose(fid);
