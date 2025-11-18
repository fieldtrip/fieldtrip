%  HBF_BINF_DIR  computes an arbitrary component of the magnetic field due to a 
%    current dipole in infinite homogeneous non-magnetic medium 
% 
%  Binf=HBF_BINF_DIR(fp,fpdir,spos,sdir)
%    fp = field points (where the field is computed), [N x 3]
%    fpdir = field orientations (which component is computed), [N x 3]
%    spos  = source positions, [M x 3]
%    spdir = dipole moments (or, if unit-norm, dipole orientations), [M x 3]
% 
%    Binf  = resulting magnetic field for each field point and dipole, [N x M]
% 
%  v160229 (c) Matti Stenroos
%