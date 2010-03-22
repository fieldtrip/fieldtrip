function write_stl(filename, pnt, tri, nrm);

% WRITE_STL writes a triangulation to an ascii *.stl file, which is a file
% format native to the stereolithography CAD software created by 3D Systems.
%
% Use as
%   write_stl(filename, pnt, tri, nrm)
% where nrm refers to the triangle normals.
%
% See also READ_STL
  
% Copyright (C) 2006, Robert Oostenveld
% 
% $Log: write_stl.m,v $
% Revision 1.1  2009/01/14 09:33:10  roboos
% moved even more files from fileio to fileio/private, see previous log entry
%
% Revision 1.1  2006/04/19 12:54:47  roboos
% new implementation, to be used in construction of solid model of monkey brain
%
  
% solid   testsphere
%   facet normal -0.13 -0.13 -0.98
%     outer loop
%       vertex 1.50000 1.50000 0.00000
%       vertex 1.50000 1.11177 0.05111
%       vertex 1.11177 1.50000 0.05111
%     endloop
%   endfacet
%   ...

fid = fopen(filename, 'wb');

ntri = size(tri,1);
fprintf(fid, 'solid unknown\n');

for i=1:ntri
  fprintf(fid, '  facet normal %f %f %f\n', nrm(i,:));
  fprintf(fid, '    outer loop\n');
  fprintf(fid, '      vertex %f %f %f\n', pnt(tri(i,1),1), pnt(tri(i,1),2), pnt(tri(i,1),3));
  fprintf(fid, '      vertex %f %f %f\n', pnt(tri(i,2),1), pnt(tri(i,2),2), pnt(tri(i,2),3));
  fprintf(fid, '      vertex %f %f %f\n', pnt(tri(i,3),1), pnt(tri(i,3),2), pnt(tri(i,3),3));
  fprintf(fid, '    endloop\n');
  fprintf(fid, '  endfacet\n');
end

fprintf(fid, 'endsolid\n');

fclose(fid);

