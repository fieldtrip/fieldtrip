function write_vtk(fn, pnt, tri)

% WRITE_VTK writes a triangulation to a VTK (Visualisation ToolKit) format file.
% Supported are triangles, tetraheders and hexaheders.
%
% Use as
%   write_vtk(filename, pnt, tri)
%
% See also READ_VTK, WRITE_PLY

% Copyright (C) 2002, Robert Oostenveld
%
% $Id$

fid = fopen_or_error(fn, 'wt');

npnt = size(pnt,1);
ntri = size(tri,1);

% write the header
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'vtk output\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');
fprintf(fid, '\n');

% write the vertex points
fprintf(fid, 'POINTS %d float\n', npnt);
fprintf(fid, '%f\t%f\t%f\n', pnt');
fprintf(fid, '\n');

if size(tri,2)==3
  % write the triangles
  fprintf(fid, 'POLYGONS %d %d\n', ntri, (3+1)*ntri);
  fprintf(fid, '3\t%d\t%d\t%d\n', (tri-1)');
  fprintf(fid, '\n');
elseif size(tri,2)==4
  % write the tetraheders
  fprintf(fid, 'POLYGONS %d %d\n', ntri, (4+1)*ntri);
  fprintf(fid, '4\t%d\t%d\t%d\t%d\n', (tri-1)');
  fprintf(fid, '\n');
elseif size(tri,2)==8
  % write the hexaheders
  fprintf(fid, 'POLYGONS %d %d\n', ntri, (8+1)*ntri);
  fprintf(fid, '8\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n', (tri-1)');
  fprintf(fid, '\n');
end

fclose(fid);
  
