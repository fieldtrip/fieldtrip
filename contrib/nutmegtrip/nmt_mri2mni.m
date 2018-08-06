function [xyz_o]=nmt_mri2mni(xyz_i,mrifullpath,doaffine)
% [xyz_o]=nmt_mri2mni(xyz_i,mrifullpath,[doaffine])
% Takes MRI coords (mm) and converts them to MNI coordinates (mm)
% using normalization info from SPM8 (or SPM12 'OldNorm')
%
% XYZ_I: Nx3 list of coords in individual's MRI coordinates (mm)
%
% mrifullpath: path to normalized MRI volume (nifti)
%
% doaffine: (optional) applies affine transform for area outside of
%           bounding box specified by SPM warping.



if size(xyz_i,2)~=3
    error('input coords must be in N x 3 matrix');
end

[mridir,mrifilebase,mrifileext] = fileparts(mrifullpath);
mrifile = [mrifilebase mrifileext];

iyimg=fullfile(mridir,['y_i' mrifile]);
% norm_mrifullpath = fullfile(mridir,['w' mrifile]);
snmat=fullfile(mridir,[mrifilebase '_sn.mat']);

if ~exist(iyimg,'file')
    disp('need to save inverse deformation fields, one moment...');
    nmt_spm_write_deformationinv(mrifullpath);
end

iyimg = [repmat(iyimg,3,1) [',1,1';',1,2';',1,3']];

% this is thanks to john's gem #5 (from now-dead link: http://www.sph.umich.edu/~nichols/JohnsGems2.html )
Viy = spm_vol(iyimg);
vx = double(nmt_transform_coord(inv(Viy(1).mat),xyz_i));  % The voxel in the deformation to sample
xyz_o = [ spm_sample_vol(Viy(1),vx(:,1),vx(:,2),vx(:,3),1) spm_sample_vol(Viy(2),vx(:,1),vx(:,2),vx(:,3),1) spm_sample_vol(Viy(3),vx(:,1),vx(:,2),vx(:,3),1)];

if ~exist('doaffine') || isempty(doaffine)
    doaffine=0;
end

if doaffine
    Vorig=spm_vol(mrifullpath);
    load(snmat);
    
    for ii=1:size(xyz_o,1)
        if (sum(xyz_o(ii,:)==0)==3 || any(isnan(xyz_o(ii,:))))
            % this will happen if voxel outside of region of Viy, and spm_sample_vol returns 'zeros' for every element
            % then we do the next best thing of applying Affine normalisation transform, rather than warped transform
            vxorig=nmt_transform_coord(inv(Vorig.mat),xyz_i(ii,:));
            xyz_o(ii,:)=nmt_transform_coord(VG.mat,nmt_transform_coord(inv(Affine),vxorig));
            disp(['warning: the ' num2str(ii) ' coordinate was out of bounds and so transformed using Affine only']);
        end
    end
end

