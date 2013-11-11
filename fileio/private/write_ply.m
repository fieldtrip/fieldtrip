function write_ply(fn, pnt, tri, format)

% WRITE_PLY writes triangles, tetraheders or hexaheders to a Stanford *.ply format file
%
% Use as
%   write_ply(filename, vertex, element)
%
% Documentation is provided on
%   http://paulbourke.net/dataformats/ply/
%   http://en.wikipedia.org/wiki/PLY_(file_format)
%
% See also READ_PLY, READ_VTK, WRITE_VTK

% Copyright (C) 2012, Robert Oostenveld
%
% $Id$

% the elements are described as
%   Mx3 list of vertices for triangles
%   Mx4 list of vertices for tetraheders
%   Mx8 list of vertices for hexaheders

if nargin<4
  format = 'ascii';
end


if strcmp(format, 'ascii')
  fid = fopen(fn, 'wt');
else
  fid = fopen(fn, 'wb');
end

if fid~=-1
  
  npnt = size(pnt,1);
  ntri = size(tri,1);
  if numel(tri)>0
    % MATLAB indexes start at 1, inside the file they start at 0
    tri = tri-1;
  end
  
  % write the header
  fprintf(fid, 'ply\n');
  if  strcmp(format, 'ascii')
    fprintf(fid, 'format ascii 1.0\n');
  else
    [c, m, e] = computer;
    if e=='L'
      fprintf(fid, 'format binary_little_endian 1.0\n');
    elseif e=='B'
      fprintf(fid, 'format binary_big_endian 1.0\n');
    end
  end
  
  fprintf(fid, 'element vertex %d\n', npnt);
  fprintf(fid, 'property float x\n');
  fprintf(fid, 'property float y\n');
  fprintf(fid, 'property float z\n');
  fprintf(fid, 'element face %d\n', ntri);
  fprintf(fid, 'property list uchar int vertex_index\n');
  fprintf(fid, 'end_header\n');
  
  if strcmp(format, 'ascii')
    % write the vertex points
    fprintf(fid, '%f\t%f\t%f\n', pnt');
    if size(tri,2)==3
      % write the triangles
      fprintf(fid, '3\t%d\t%d\t%d\n', tri');
    elseif size(tri,2)==4
      % write the tetraheders
      fprintf(fid, '4\t%d\t%d\t%d\t%d\n', tri');
    elseif size(tri,2)==8
      % write the hexaheders
      fprintf(fid, '8\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n', tri');
    end
  else
    fwrite(fid, pnt', 'single');      % float
    num = size(tri,2);
    for i=1:size(tri,1)
      fwrite(fid, num, 'uint8');      % uchar
      fwrite(fid, tri(i,:), 'int32'); % int
    end
  end
  
  fclose(fid);
  
else
  error('unable to open file');
end
