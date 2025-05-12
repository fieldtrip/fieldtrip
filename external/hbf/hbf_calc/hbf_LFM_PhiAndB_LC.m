function [LFMphi,LFMb] = hbf_LFM_PhiAndB_LC(bmeshes,coils,Tphi,TBvol,spos,sdir,averef)
% HBF_LFM_PHIANDB_LC builds electric and magnetic lead field matrices
%
% LFM=HBF_LFM_PHIANDB_LC(bmeshes,coils,Tphi,TBvol,spos,sdir,averef)
% LFM=HBF_LFM_PHIANDB_LC(bmeshes,coils,Tphi,TBvol,spos,averef)
%
% You can compute potential due to any set of directed dipoles by
% giving the dipole moments (with amplitude) in the 'sdir' argument. If you
% omit 'sdir', LFM will be computed using xyz unit dipole triplets
%
%   bmeshes: BEM geometry, cell array of hbf mesh structs
%   coils:  coil description, hbf struct
%   Tphi:   BEM transfer matrix built with hbf tools
%   TBvol:  TBvol matrix built with the hbf BEM solver
%   spos:   source positions, [M x 3]
%   sdir:   source orientations (unit-length) or moments, [M x 3]
%   averef: use average reference: 1 for yes (default), 2 for no.
%
% If sdir is not given,
%   LFMphi: electric lead field matrix, [Number of electrodes x 3M]
%       [tphi_1x tphi_1y tphi_1z ... tphi_Mx tphi_My tphi_Mz]
%   LFMb: magnetic lead field matrix, [Number of coils x 3M]
%       [tb_1x tb_1y tb_1z ... tb_Mx tb_My tb_Mz]
%
% If sdir is given,
%   LFMphi: electric lead field matrix, [Number of electrodes x M]
%       [tphi_1 ... tphi_M]
%   LFMb: magnetic lead field matrix, [Number of coils x M]
%       [tb_1 ... tb_M]
%
%
% This function assumes average reference for potential by default. If
% 'flag_averef' is selected 0, the potential is computed againts the
% reference chosen when building the BEM matrix (for ISA approach, the
% average over the ISA surface; for non-ISA BEM, the outermost surface).
%
% This is a wrapper for HBF_LFM_PHI_LC_DIR and HBF_LFM_PHI_LC_XYZ
%
% v200928 (c) Matti Stenroos

% parse for the use of average reference
averef_flag = 1;
if nargin == 7
    averef_flag = averef;
elseif nargin == 6 && numel(sdir) == 1
    averef_flag = sdir;
end

if nargin == 7 || (nargin == 6 && size(sdir,2) == 3)
    fprintf('Computing oriented lead fields for Phi and B...');
    [LFMphi,LFMb] = hbf_LFM_PhiAndB_LC_dir(bmeshes,coils,Tphi,TBvol,spos,sdir,averef_flag);
else
    fprintf('Computing xyz lead fields for Phi and B...');
    [LFMphi,LFMb] = hbf_LFM_PhiAndB_LC_xyz(bmeshes,coils,Tphi,TBvol,spos,averef_flag);
end
fprintf('OK.\n');