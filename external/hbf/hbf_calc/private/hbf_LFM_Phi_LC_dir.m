%  HBF_LFM_PHI_LC_DIR builds electric lead field matrix based on directed 
%  current dipoles.
% 
%  LFM=HBF_LFM_PHI_LC_DIR(meshes,Tphi,spos,sdir,flag_averef)
%    meshes: BEM geometry, cell array of hbf structs
%    Tphi:   Tphi matrix built with the hbf BEM solver
%    spos:   source positions, [M x 3]
%    sdir:   source orientations (unit-length), [M x 3]
%    flag_averef (optional, default value 1): give 0, if you do not want to
%            use average reference
% 
%    LFM:   lead field matrix, [Number of electrodes x M]
%        [t_1 ... t_M]
% 
%  You can also compute phi due to any set of directed dipoles by 
%  giving the dipole moments (with amplitude) in the 'sdir' argument.
% 
%  This function assumes average reference by default. If 'flag_averef' is
%  selected 0, the potential is computed againts the reference chosen when
%  building the BEM matrix.
% 
%  v160229 (c) Matti Stenroos
%