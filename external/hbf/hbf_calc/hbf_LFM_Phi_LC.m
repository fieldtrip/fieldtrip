function LFM = hbf_LFM_Phi_LC(bmeshes,Tphi,spos,sdir,averef)
% HBF_LFM_Phi_LC builds the electric lead field matrix
%
% LFM=HBF_LFM_PHI_LC(bmeshes,Tphi,spos,sdir,averef)
% LFM=HBF_LFM_PHI_LC(bmeshes,Tphi,spos,averef)
%
% You can compute potential due to any set of directed dipoles by 
% giving the dipole moments (with amplitude) in the 'sdir' argument. If you
% omit 'sdir', LFM will be computed using xyz unit dipole triplets
%
%   bmeshes: BEM geometry, cell array of hbf mesh structs
%   coils:  coil description, hbf struct
%   Tphi:   BEM transfer matrix built with hbf tools
%   spos:   source positions, [M x 3]
%   sdir:   source orientations (unit-length) or moments, [M x 3]
%   averef: use average reference: 1 for yes (default), 2 for no.
%
%   LFM:   lead field matrix, [Number of electrodes x M]
%           [t_1 ... t_M]
%           or, if 'sdir' omitted, [Number of electrodes x 3M]
%           [t_1x t_1y t1_z ... t_Mx t_My t_Mz]
%
% If you set averef = 0, the lead field matrix is computed against the
% reference set when building the BEM model (for ISA approach, average
% over the ISA surface; for non-ISA BEM, the outermost surface).
%
% This is a wrapper for HBF_LFM_PHI_LC_DIR and HBF_LFM_PHI_LC_XYZ 
%
% v200928 (c) Matti Stenroos

% parse for the use of average reference
averef_flag = 1;
if nargin == 5
    averef_flag = averef;
elseif nargin == 4 && numel(sdir) == 1
    averef_flag = sdir;
end

if nargin == 5 || (nargin == 4 && size(sdir,2) == 3)
    fprintf('Computing oriented lead fields for Phi...');
    LFM = hbf_LFM_Phi_LC_dir(bmeshes,Tphi,spos,sdir,averef_flag);
else
    fprintf('Computing xyz lead fields for Phi...');
    LFM = hbf_LFM_Phi_LC_xyz(bmeshes,Tphi,spos,averef_flag);
end
fprintf('OK.\n');