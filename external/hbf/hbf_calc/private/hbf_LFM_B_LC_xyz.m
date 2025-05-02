%  HBF_LFM_B_LC_XYZ builds magnetic lead field matrix based on xyz-oriented
%    unit-current dipoles
% 
%  LFMm=HBF_LFM_B_LC_XYZ(meshes,coils,TB,spos)
%    meshes: BEM geometry, cell array of hbf structs
%    coils:  coil description, hbf struct
%    TB:     TB matrix built with the hbf BEM solver
%    spos:   source positions, [M x 3]
% 
%    LFMm:   lead field matrix, [Number of coils x 3M]
%        [t_1x t_1y t1_z ... t_Mx t_My t_Mz]
% 
%  v160229 (c) Matti Stenroos
%