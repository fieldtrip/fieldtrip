function write_vtk(fn, pnt, dhk)

% WRITE_VTK writes a triangulation to a VTK (Visualisation ToolKit) format file.
% Supported are triangles, tetraheders and hexaheders.
%
% Use as
%   write_vtk(filename, pnt, dhk)
%
% See also READ_VTK, WRITE_PLY

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
  
  if size(dhk,2)==3
    % write the triangles
    fprintf(fid, 'POLYGONS %d %d\n', ndhk, (3+1)*ndhk);
    fprintf(fid, '3\t%d\t%d\t%d\n', (dhk-1)');
    fprintf(fid, '\n');
  elseif size(dhk,2)==4
    % write the tetraheders
    fprintf(fid, 'POLYGONS %d %d\n', ndhk, (4+1)*ndhk);
    fprintf(fid, '4\t%d\t%d\t%d\t%d\n', (dhk-1)');
    fprintf(fid, '\n');
  elseif size(dhk,2)==8
    % write the hexaheders
    fprintf(fid, 'POLYGONS %d %d\n', ndhk, (8+1)*ndhk);
    fprintf(fid, '8\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n', (dhk-1)');
    fprintf(fid, '\n');
  end
  
  fclose(fid);
  
else
  error('unable to open file');
end
