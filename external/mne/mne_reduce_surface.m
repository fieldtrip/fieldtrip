function [verts, faces] = mne_reduce_surface(surfin,desired_ntri,surfout)
%
%  [verts, faces] = mne_reduce_surface(surfin,desired_ntri,surfout)
%
%  verts       - Vertex coordinates in meters
%  faces       - The triangulation information
%
%  surfin      - Name of a surface file to read
%  desired_nri - Desired number of triangles after reduction
%  surfout     - Name of a surface file to hold the reduce surface (optional)
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.2  2009/02/03 10:26:13  msh
%   Added reading of 'new' quadrangle files
%
%   Revision 1.1  2009/01/28 20:43:52  msh
%   Added mne_reduce_surface
%
%
me='MNE:mne_reduce_surface';

if nargin < 2
   error(me,'Incorrect number of input arguments');
end

fprintf(1,'Reading %s...\n',surfin);
[surf.vertices,surf.faces] = mne_read_surface(surfin);

ratio = desired_ntri/size(surf.faces,1);

fprintf(1,'Reducing the number of triangles...');
surf2 = reducepatch(surf,ratio);
fprintf(1,'[done]\n');

fprintf(1,'After reduction: %d vertices and %d triangles\n',size(surf2.vertices,1),size(surf2.faces,1));

if nargin == 3
   mne_write_surface(surfout,surf2.vertices,surf2.faces);
end

verts = surf2.vertices;
faces = surf2.faces;

return
