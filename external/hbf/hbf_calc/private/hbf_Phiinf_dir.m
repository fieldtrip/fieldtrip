%  HBF_PHIINF_DIR  computes potential due to a current dipole in infinite
%    homogeneous conductor that has unit conductivity
% 
%  phiinf=HBF_PHIINF_DIR(fp,spos,smom)
%    fp = field points (where the field is computed), [N x 3]
%    spos  = source positions, [M x 3]
%    smom = dipole moments, [M x 3]
% 
%    phiinf = resulting potential for each field point and dipole, [N x M]
% 
%  v160229 (c) Matti Stenroos
%