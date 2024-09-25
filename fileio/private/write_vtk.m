function write_vtk(fn, pos, tri, val)

% WRITE_VTK writes a mesh to a VTK (Visualisation ToolKit) format file.
% Supported are triangles, tetraheders and hexaheders.
%
% Use as
%   write_vtk(filename, pos, tri, val)
%   write_vtk(filename, pos, tet, val)
%   write_vtk(filename, pos, hex, val)
% where pos describes the vertex positions and tri/tet/hex describe the connectivity
% of the surface or volume elements. 
% 
% The optional val argument can be used to write scalar or vector values for
% each vertex or element.
%
% See also READ_VTK, WRITE_PLY

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

if nargin<4
  val = [];
end

npos = size(pos,1);
ntri = size(tri,1);   % this can be triangles, tetraheders or hexaheders

if ~isempty(val)
  assert(size(val,1)==npos || size(val,1)==ntri, 'the values are inconsistent with the mesh dimensions');
end

% see https://kitware.github.io/vtk-examples/site/Cxx/GeometricObjects/LinearCellDemo/
VTK_VERTEX           = 1;
VTK_POLY_VERTEX      = 2;
VTK_LINE             = 3;
VTK_POLY_LINE        = 4;
VTK_TRIANGLE         = 5;
VTK_TRIANGLE_STRIP   = 6;
VTK_POLYGON          = 7;
VTK_PIXEL            = 8;
VTK_QUAD             = 9;
VTK_TETRA            = 10;
VTK_VOXEL            = 11;
VTK_HEXAHEDRON       = 12;
VTK_WEDGE            = 13;
VTK_PYRAMID          = 14;
VTK_PENTAGONAL_PRISM = 15;
VTK_HEXAGONAL_PRISM  = 16;

fid = fopen_or_error(fn, 'wt');

% write the header
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'vtk output\n');
fprintf(fid, 'ASCII\n');

if size(tri,2)==3
  if isempty(val)
    fprintf(fid, 'DATASET POLYDATA\n');
    fprintf(fid, '\n');
    % write the vertices
    fprintf(fid, 'POINTS %d float\n', npos);
    fprintf(fid, '%f\t%f\t%f\n', pos');
    fprintf(fid, '\n');
    % write the triangles
    fprintf(fid, 'POLYGONS %d %d\n', ntri, (3+1)*ntri);
    fprintf(fid, '3\t%d\t%d\t%d\n', (tri-1)');
    fprintf(fid, '\n');
  else
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, '\n');
    % write the vertices
    fprintf(fid, 'POINTS %d float\n', npos);
    fprintf(fid, '%f\t%f\t%f\n', pos');
    fprintf(fid, '\n');
    % write the triangles
    fprintf(fid, 'CELLS %d %d\n', ntri, (3+1)*ntri);
    fprintf(fid, '3\t%d\t%d\t%d\n', (tri-1)');
    fprintf(fid, '\n');
    fprintf(fid, 'CELL_TYPES %d\n', ntri);
    fprintf(fid, '%d\n', repmat(VTK_TRIANGLE, ntri, 1));
    fprintf(fid, '\n');
  end

elseif size(tri,2)==4
  fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
  fprintf(fid, '\n');
  % write the vertices
  fprintf(fid, 'POINTS %d float\n', npos);
  fprintf(fid, '%f\t%f\t%f\n', pos');
  fprintf(fid, '\n');
  % write the tetraheders
  fprintf(fid, 'CELLS %d %d\n', ntri, (4+1)*ntri);
  fprintf(fid, '4\t%d\t%d\t%d\t%d\n', (tri-1)');
  fprintf(fid, '\n');
  fprintf(fid, 'CELL_TYPES %d\n', ntri);
  fprintf(fid, '%d\n', repmat(VTK_TETRA, ntri, 1));
  fprintf(fid, '\n');

elseif size(tri,2)==8
  fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
  fprintf(fid, '\n');
  % write the vertices
  fprintf(fid, 'POINTS %d float\n', npos);
  fprintf(fid, '%f\t%f\t%f\n', pos');
  fprintf(fid, '\n');
  % write the hexaheders
  fprintf(fid, 'CELLS %d %d\n', ntri, (8+1)*ntri);
  fprintf(fid, '8\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n', (tri-1)');
  fprintf(fid, '\n');
  fprintf(fid, 'CELL_TYPES %d\n', ntri);
  fprintf(fid, '%d\n', repmat(VTK_HEXAHEDRON, ntri, 1));
  fprintf(fid, '\n');

else
  error('cannot write this type of data')
end

if size(val,1)==size(pos,1)
  % write scalar values corresponding with the vertices
  fprintf(fid, 'POINT_DATA %d\n', ntri);
  fprintf(fid, 'SCALARS data float %d\n', size(val,2));
  fprintf(fid, 'LOOKUP_TABLE default\n');
  % the format depends on the number of items per vertex
  fmt = '%f\t';
  fmt = repmat(fmt, 1, size(val,2));
  fmt(end) = 'n';
  fprintf(fid, fmt, val');
elseif size(val,1)==size(tri,1)
  % write scalar values corresponding with the hexaheders
  fprintf(fid, 'CELL_DATA %d\n', ntri);
  fprintf(fid, 'SCALARS data float %d\n', size(val,2));
  fprintf(fid, 'LOOKUP_TABLE default\n');
  % the format depends on the number of items per vertex
  fmt = '%f\t';
  fmt = repmat(fmt, 1, size(val,2));
  fmt(end) = 'n';
  fprintf(fid, fmt, val');
end

fclose(fid);
