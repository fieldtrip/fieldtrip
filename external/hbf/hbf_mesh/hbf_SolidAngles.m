% HBF_SOLIDANGLES calculates solid angles spanned by triangles of a mesh in
%    one point or a pointset
% 
% omega=HBF_SOLIDANGLES(mesh,fieldpoints)
% omega=HBF_SOLIDANGLES(elements,vertices,fieldpoints)
%    mesh:       hbf mesh struct
%    elements:   triangle description, [N(triangles) x 3]
%    vertices:   mesh vertices, [N(vertices) x 3]
%    fieldpoint: point, in which the angles are evaluated, [N(fieldpoints) x 3]
% 
%    omegas: solid angles, [N(fieldpoints) x N(triangles]
%  Eq. 8 of vanOosterom & Strackee, IEEE TBME 1983
% 
%  This function assumes, by tradition, that the mesh has CW orientation.
%  When used with CCW meshes (like meshes in Helsinki BEM framework), the
%  angles thus have wrong sign: a point inside a closed mesh sees omega of
%  -4pi instead of correct 4pi.
% 
%  v160419 (c) Matti Stenroos
%