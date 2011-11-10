function write_tri(fn, pnt, dhk)

% WRITE_VTK writes a triangulation to a VTK (Visualisation ToolKit) format file
%
% write_vtk(filename, pnt, dhk)
%
% See also READ_VTK, READ_TRI, WRITEE_TRI

% Copyright (C) 2002, Robert Oostenveld
%
% $Id$
fid = fopen(fn, 'wt');
if fid~=-1

  npnt = size(pnt,1);
  ndhk = size(dhk,1);

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

  % write the triangles
  fprintf(fid, 'POLYGONS %d %d\n', ndhk, 4*ndhk);
  fprintf(fid, '3\t%d\t%d\t%d\n', (dhk-1)');
  fprintf(fid, '\n');

  fclose(fid);

else
  error('unable to open file');
end


