function LFM=hbf_LFM_B_LC(bmeshes,coils,TB,spos,sdir)
% HBF_LFM_B_LC builds the magnetic lead field matrix
%
% LFM=HBF_LFM_B_LC(bmeshes,coils,TB,spos,sdir)
% LFM=HBF_LFM_B_LC(bmeshes,coils,TB,spos)
%
% You can also compute magnetic field due to any set of directed dipoles by 
% giving the dipole moments (with amplitude) in the 'sdir' argument. If you
% omit 'sdir', LFM will be computed using xyz unit dipole triplets
%
%   bmeshes: BEM geometry, cell array of hbf mesh structs
%   coils:  coil description, hbf struct
%   TBvol:  TBvol matrix built with the hbf BEM solver
%   spos:   source positions, [M x 3]
%   sdir:   source orientations (unit-length), [M x 3]; optional
%
%   LFM:   lead field matrix, [Number of coils (field points) x M]
%           [t_1 ... t_M]
%           or, if 'sdir' omitted, [Number of coils x 3M]
%           [t_1x t_1y t1_z ... t_Mx t_My t_Mz]
%
% This is a wrapper for HBF_LFM_B_LC_DIR and HBF_LFM_B_LC_XYZ 
%
% v200924 (c) Matti Stenroos
if nargin==5
    fprintf('Computing oriented lead fields for B...');
    LFM=hbf_LFM_B_LC_dir(bmeshes,coils,TB,spos,sdir);
else
    fprintf('Computing xyz lead fields for B...');
    LFM=hbf_LFM_B_LC_xyz(bmeshes,coils,TB,spos);
end
fprintf('OK.\n');