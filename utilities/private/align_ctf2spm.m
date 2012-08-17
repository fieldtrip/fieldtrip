function [mri] = align_ctf2spm(mri)

% ALIGN_CTF2SPM performs an approximate alignment of the anatomical volume
% from CTF towards SPM coordinates. Only the homogeneous transformation matrix
% is modified and the coordsys-field is updated.

%--------------------------------------------------------------------------
% do a first round of approximate coregistration using the predefined voxel
% locations of the fiducials and landmarks in the template T1 image from
% SPM

spmvox2spmhead = [
     2     0     0   -92
     0     2     0  -128
     0     0     2   -74
     0     0     0     1
];

% these are the voxel indices of some points in the SPM canonical T1 
spmvox_Ac           = [46 64  37  1]';     % the anterior commissure
spmvox_Ori          = [46 48  10  1]';     % approximately between the ears in T1.mnc
spmvox_Nas          = [46 106 13  1]';     % approximately the nasion in T1.mnc
spmvox_Lpa_canal    = [ 5 48  10  1]';     % Left ear canal
spmvox_Rpa_canal    = [87 48  10  1]';     % Right ear canal

spmhead_Ac           = spmvox2spmhead * spmvox_Ac  ;
spmhead_Ori          = spmvox2spmhead * spmvox_Ori ;
spmhead_Nas          = spmvox2spmhead * spmvox_Nas ;
spmhead_Lpa_canal    = spmvox2spmhead * spmvox_Lpa_canal ;
spmhead_Rpa_canal    = spmvox2spmhead * spmvox_Rpa_canal ;

ctfvox2ctfhead  = mri.transform;
spmhead2ctfhead = headcoordinates(spmhead_Nas(1:3), spmhead_Lpa_canal(1:3), spmhead_Rpa_canal(1:3), 'ctf');

%ctfvox2spmhead =  inv(spmhead2ctfhead) *  ctfvox2ctfhead;
ctfvox2spmhead =  spmhead2ctfhead \ ctfvox2ctfhead;

% change the transformation matrix, such that it returns approximate SPM head coordinates
mri.transform = ctfvox2spmhead;
mri.coordsys  = 'spm';

% do a second round of affine registration (rigid body) to get improved
% alignment with spm coordinate system. this is needed because there may be
% different conventions defining LPA and RPA.
switch spm('ver')
  case 'SPM8'
    template = fullfile(spm('Dir'),'templates','T1.nii');
  case 'SPM2'
    template = fullfile(spm('Dir'),'templates','T1.mnc');
  otherwise
    error('unsupported spm-version');    
end
mri2 = ft_read_mri(template);

tname1 = [tempname, '.img'];
tname2 = [tempname, '.img'];
V1 = ft_write_mri(tname1, mri.anatomy,  'transform', mri.transform,  'spmversion', spm('ver'));
V2 = ft_write_mri(tname2, mri2.anatomy, 'transform', mri2.transform, 'spmversion', spm('ver'));

flags.regtype = 'rigid';
[M, scale]    = spm_affreg(V1,V2,flags);

% some juggling around with the transformation matrices
ctfvox2spmhead2  = M \ V1.mat;
spmhead2ctfhead2 = ctfvox2ctfhead / ctfvox2spmhead2;

% append the transformation matrices to the output
mri.transform     = ctfvox2spmhead2;
mri.vox2headOrig  = ctfvox2ctfhead;
mri.vox2head      = ctfvox2spmhead2;
mri.head2headOrig = spmhead2ctfhead2;

% delete the temporary files
delete(tname1);
delete(tname2);
