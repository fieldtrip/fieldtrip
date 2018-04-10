function [xyz_o_mm,xyz_o_vx]=nmt_mni2mri(xyz_i,mrifullpath)
% [xyz_o_mm,xyz_o_vx]=nmt_mri2mni(xyz_i,mrifullpath,[doaffine])
% Takes MNI coordinates (mm) and converts to original MRI coordinates (mm
% and/or voxel) using normalization info from SPM8 (or SPM12 'OldNorm')
%
% XYZ_I: Nx3 list of coords in individual's MRI coordinates (mm)
% mrifullpath: path to normalized MRI volume (nifti)
%
% [wrapper for spm_get_orig_coord from SPM8/SPM12]

if size(xyz_i,2)~=3
    error('input coords must be in N x 3 matrix');
end

[mridir,mrifilebase,mrifileext] = fileparts(mrifullpath);
mrifile = [mrifilebase mrifileext];

snmat=fullfile(mridir,[mrifilebase '_sn.mat']);

xyz_o_mm = spm_get_orig_coord(xyz_i,snmat);
xyz_o_vx = spm_get_orig_coord(xyz_i,snmat,mrifullpath);