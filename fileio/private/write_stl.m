function write_stl(filename, pnt, tri, nrm)

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
  
% solid   testsphere
%   facet normal -0.13 -0.13 -0.98
%     outer loop
%       vertex 1.50000 1.50000 0.00000
%       vertex 1.50000 1.11177 0.05111
%       vertex 1.11177 1.50000 0.05111
%     endloop
%   endfacet
%   ...

fid = fopen_or_error(filename, 'wb');

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

