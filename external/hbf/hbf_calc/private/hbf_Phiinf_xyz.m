%  HBF_PHIINF_XYZ  computes potential due to a xyz-oriented unit-current-dipole 
%    triplet in infinite homogeneous conductor that has unit conductivity
% 
%  phiinf=HBF_PHIINF_XYZ(fp,spos)
%    fp = field points (where the field is computed), [N x 3]
%    spos  = source positions, [M x 3]
% 
%    phiinf = resulting potential for each field point and dipole, [N x 3M]
%      [phiinf_1x phiinf_1y phiinf_1z ... phiinf_Mx phiinf_My phiinf_Mz]
% 
%  v160229 (c) Matti Stenroos
%