function [pnt, tri, srf] = read_bv_srf(filename);

% READ_BV_SRF reads a triangulated surface from a BrainVoyager *.srf file
%
% Use as
%   [pnt, tri] = read_bv_srf(filename) or
%   [pnt, tri, srf] = read_bv_srf(filename)

% Copyright (C) 2005, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% This documentation originates from
%   http://www.brainvoyager.com/BV2000OnlineHelp/BrainVoyagerWebHelp/mergedProjects/FileFormats/BrainVoyager_File_Formats.htm
% 
% BYTES DATA TYPE DEFAULT DESCRIPTION
% 4 float 3 version number
% 4 int 0 reserved, must be '0'
% 4 int  NrOfVertices (number of vertices)
% 4 int  NrOfTriangles (number of triangles)
% 4 float 128.0 MeshCenterX
% 4 float 128.0 MeshCenterY
% 4 float 128.0 MeshCenterZ
% NrOfVertices*4 float  VertexX, sequence of X coordinates of all vertices
% NrOfVertices*4 float  VertexY, sequence of Y coordinates of all vertices
% NrOfVertices*4 float  VertexZ, sequence of Z coordinates of all vertices
% NrOfVertices*4 float  NormalX, sequence of X components of all vertex normals
% NrOfVertices*4 float  NormalY, sequence of Y components of all vertex normals
% NrOfVertices*4 float  NormalZ, sequence of Z components of all vertex normals
% 4 float 0.322 R component of convex curvature color (range: 0.0 - 1-0)
% 4 float 0.733 G component of convex curvature color (range: 0.0 - 1-0)
% 4 float 0.980 B component of convex curvature color (range: 0.0 - 1-0)
% 4 float 0.500 Alpha component of convex curvature color (range: 0.0 - 1-0)
% 4 float 0.100 R component of concave curvature color (range: 0.0 - 1-0)
% 4 float 0.240 G component of concave curvature color (range: 0.0 - 1-0)
% 4 float 0.320 B component of concave curvature color (range: 0.0 - 1-0)
% 4 float 0.500 Alpha component of concave curvature color (range: 0.0 - 1-0)
% NrOfVertices*4 int  MeshColor, sequence of color indices of all vertices (*1)
% 4 int  N, number of nearest neighbors of vertex 1
% 4 int  Nearest neighbor 1 of vertex 1
% : :  :
% 4 int  Nearest neighbor N of vertex 1
% :: ::  ::
% 4 int  N, number of nearest neighbors of vertex 'NrOfVertices'
% 4 int  Nearest neighbor 1 of vertex 'NrOfVertices'
% : :  :
% 4 int  Nearest neighbor N of vertex 'NrOfVertices'
% NrOfTriangels*3*4 int  Sequence of three indices to constituting vertices of each triangle
% 4 int  NrOfTriangleStripElements
% NrOfStripElements*4 int  Sequence of strip elements (if NrOfStripElements > 0)


fid = fopen(filename, 'rb', 'ieee-le');

srf.version_number                                        = fread(fid, 1, 'float');
srf.reserved                                              = fread(fid, 1, 'int'  );
srf.NrOfVertices                                          = fread(fid, 1, 'int'  );
srf.NrOfTriangles                                         = fread(fid, 1, 'int'  );
NrOfVertices  = srf.NrOfVertices;
NrOfTriangles = srf.NrOfTriangles;
srf.MeshCenterX                                           = fread(fid, 1, 'float');
srf.MeshCenterY                                           = fread(fid, 1, 'float');
srf.MeshCenterZ                                           = fread(fid, 1, 'float');
srf.X_coordinates                                         = fread(fid, NrOfVertices, 'float');
srf.Y_coordinates                                         = fread(fid, NrOfVertices, 'float');
srf.Z_coordinates                                         = fread(fid, NrOfVertices, 'float');
srf.X_components                                          = fread(fid, NrOfVertices, 'float');
srf.Y_components                                          = fread(fid, NrOfVertices, 'float');
srf.Z_components                                          = fread(fid, NrOfVertices, 'float');
srf.R_component_of_convex_curvature_color                 = fread(fid, 1, 'float');
srf.G_component_of_convex_curvature_color                 = fread(fid, 1, 'float');
srf.B_component_of_convex_curvature_color                 = fread(fid, 1, 'float');
srf.Alpha_component_of_convex_curvature_color             = fread(fid, 1, 'float');
srf.R_component_of_concave_curvature_color                = fread(fid, 1, 'float');
srf.G_component_of_concave_curvature_color                = fread(fid, 1, 'float');
srf.B_component_of_concave_curvature_color                = fread(fid, 1, 'float');
srf.Alpha_component_of_concave_curvature_color            = fread(fid, 1, 'float');
srf.MeshColor                                             = fread(fid, NrOfVertices, 'int'  );
for i=1:NrOfVertices
  number           = fread(fid, 1, 'int'  );
  srf.neighbour{i} = fread(fid, number, 'int'  );
end
srf.Triangles                                             = fread(fid, [3 NrOfTriangles], 'int'  );
srf.NrOfTriangleStripElements                             = fread(fid, 1, 'int'  );
srf.sequence_of_strip_elements                            = fread(fid, srf.NrOfTriangleStripElements, 'int'  );

fclose(fid);

pnt = [srf.X_coordinates(:) srf.Y_coordinates(:) srf.Z_coordinates(:)];
tri = srf.Triangles' + 1;

