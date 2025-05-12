%  HBF_LFM_PHIANDB_LC_DIR builds electric and magnetic lead field matrices
%    based on directed current dipoles.
% 
%  [LFMphi,LFMm]=HBF_LFM_PHIANDB_LC_DIR(meshes,coils,Tphi,TB,spos,sdir,flag_averef)
%    meshes: BEM geometry, cell array of hbf structs
%    coils:  coil description, hbf struct
%    TB:     TB matrix built with the hbf BEM solver
%    Tphi:   Tphi matrix built with the hbf BEM solver
%    spos:   source positions, [M x 3]
%    sdir:   source orientations (unit-length) or dipole moments, [M x 3]
%    flag_averef (optional, default value 1): give 0, if you do not want to
%            use average reference
% 
%    LFMphi: electric lead field matrix, [Number of electrodes x M]
%        [tphi_1 ... tphi_M]
%    LFMb: magnetic lead field matrix, [Number of coils x M]
%        [tb_1 ... tb_M]
% 
% 
%  You can compute phi and B due to any set of directed dipoles by
%  giving the dipole moments (with amplitude) in the 'sdir' argument.
% 
%  This function assumes average reference for potential by default. If
%  'flag_averef' is selected 0, the potential is computed againts the
%  reference chosen when building the BEM matrix (for ISA approach, the
%  average over the ISA surface; for non-ISA BEM, the outermost surface).
% 
%  v160229 (c) Matti Stenroos
%