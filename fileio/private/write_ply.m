function write_ply(fn, pnt, dhk)

% WRITE_PLY writes triangles, tetraheders or hexaheders to a Stanford *.ply format file
%
% Use as
%   write_ply(filename, vertex, element)
%
% Documentation is provided on
%   http://paulbourke.net/dataformats/ply/
%   http://en.wikipedia.org/wiki/PLY_(file_format)
%
% See also WRITE_VTK, READ_VTK

% Copyright (C) 2012, Robert Oostenveld
%
% $Id$

fid = fopen(fn, 'wt');
if fid~=-1
  
  npnt = size(pnt,1);
  ndhk = size(dhk,1);
  
  % write the header
  fprintf(fid, 'ply\n');
  fprintf(fid, 'format ascii 1.0\n');
  fprintf(fid, 'element vertex %d\n', npnt);
  fprintf(fid, 'property float x\n');
  fprintf(fid, 'property float y\n');
  fprintf(fid, 'property float z\n');
  fprintf(fid, 'element face %d\n', ndhk);
  fprintf(fid, 'property list uchar int vertex_index\n');
  fprintf(fid, 'end_header\n');
  
  % write the vertex points
  fprintf(fid, '%f\t%f\t%f\n', pnt');
  
  if size(dhk,2)==3
    % write the triangles
    fprintf(fid, '3\t%d\t%d\t%d\n', (dhk-1)');
  elseif size(dhk,2)==4
    % write the tetraheders
    fprintf(fid, '4\t%d\t%d\t%d\t%d\n', (dhk-1)');
  elseif size(dhk,2)==8
    % write the hexaheders
    fprintf(fid, '8\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n', (dhk-1)');
  end
  
  fclose(fid);
  
else
  error('unable to open file');
end

