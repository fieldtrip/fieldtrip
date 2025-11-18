%  HBF_BINF_XYZ  computes an arbitrary component of the magnetic field due to 
%    xyz-oriented unit-current-dipole triplet in infinite homogeneous non-magnetic medium 
% 
%  Binf=HBF_BINF_XYZ(fp,fpdir,spos)
%    fp = field points (where the field is computed), [N x 3]
%    fpdir = field orientations (which component is computed), [N x 3]
%    spos  = source positions, [M x 3]
% 
%    Binf  = resulting magnetic field for each field point and dipole, [N x 3M]
%        [B_1x B_1y B_1z ... B_Mx B_My B_Mz]
% 
%  v160229 (c) Matti Stenroos
%