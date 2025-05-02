% HBF_PROJECTPOINTSTOMESH ... projects a set of points onto a mesh.
% 
%  [pproj,loctype,locinfo,projdist]=HBF_PROJECTPOINTSTOMESH(mesh,pointset)
% 
%  This function projects the points of a pointset to the nearest positions 
%  on the mesh.
% 
%    mesh: triangle mesh, hbf struct
%    pointset: set of points to be projected, [N x 3]
%    pproj: the projected points, [N x 3]
%    loctype: location information of the projection point, [N x 1],
%        1 = in a triangle, 2 = on triangle-side, 3 = in node
%    locinfo: further location information depending on loctype, [N x 2], 
%        triangle index, two node indices, or one node index
%    projdist: distance between the projected point and original point, [N x 1] 
%  
%    The algorithm tests the distance of each point to all nodes, sides and
%    triangles of the mesh, and thus always finds the nearest points. It is
%    thus as an algorithm rather slow, but the implementation is quite well
%    optimized.
% 
%  v200331 (c) Matti Stenroos, matti.stenroos@aalto.fi
%