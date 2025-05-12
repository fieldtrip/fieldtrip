%  HBF_D_LC makes a LC double-layer matrix between two meshes
% 
%  D=HBF_D_LC(fieldmesh,sourcemesh,verbose)
% 
%    fieldmesh:  mesh, where double-layer potential is evaluated, hbf struct
%    sourcemesh: mesh, where the double layer is spanned, hbf struct
%    verbose:    optional intermediate output. 1 for yes, 0 for no (default)
% 
%    D: Double-layer matrix, [N(field vertices) x N(source vertices)]
% 
%  v190403 (c) Matti Stenroos
%