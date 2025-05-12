%  HBF_LFM_PHIANDB_LC_XYZ builds electric and magnetic lead field matrices
%    based on xyz-oriented unit-current dipoles.
% 
%  [LFMphi,LFMm]=HBF_LFM_PHIANDB_LC_XYZ(meshes,coils,Tphi,TB,spos,flag_averef)
%    meshes: BEM geometry, cell array of hbf structs
%    coils:  coil description, hbf struct
%    Tphi:   Tphi matrix built with the hbf BEM solver
%    TBvol:  TBvol matrix built with the hbf BEM solver
%    spos:   source positions, [M x 3]
%    flag_averef (optional, default value 1): give 0, if you do not want to
%            use average reference in LFMphi
% 
%    LFMphi: electric lead field matrix, [Number of electrodes x 3M]
%        [tphi_1x tphi_1y tphi_1z ... tphi_Mx tphi_My tphi_Mz]
%    LFMb: magnetic lead field matrix, [Number of coils x 3M]
%        [tb_1x tb_1y tb_1z ... tb_Mx tb_My tb_Mz]
% 
%  This function assumes average reference by default. If 'flag_averef' is
%  selected 0, the potential is computed againts the reference chosen when
%  building the BEM matrix (for ISA approach, the average over the ISA
%  surface; for non-ISA BEM, the outermost surface).
% 
%  v160303 (c) Matti Stenroos
%