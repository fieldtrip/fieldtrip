function [mri] = align_itab2spm(mri, opt)

% ALIGN_ITAB2SPM performs an approximate alignment of the anatomical volume
% from ITAB towards SPM coordinates. Only the homogenous transformation matrix
% is modified and the coordsys-field is updated.
%
% Use as
%   mri = align_itab2spm(mri)
%   mri = align_itab2spm(mri, opt)
%
% Where mri is a FieldTrip MRI-structure, and opt an optional argument
% specifying how the registration is done.
%   opt = 0: only an approximate coregistration
%   opt = 1: an approximate coregistration, followed by spm_affreg
%   opt = 2: (default): an approximate coregistration, followed by spm_normalise

if nargin<2
  opt = 2;
end

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

itabvox2itabhead  = mri.transform;
spmhead2itabhead  = ft_headcoordinates(spmhead_Nas(1:3), spmhead_Lpa_canal(1:3), spmhead_Rpa_canal(1:3), 'itab');

%itabvox2spmhead =  inv(spmhead2itabhead) * itabvox2itabhead;
itabvox2spmhead  =  spmhead2itabhead \ itabvox2itabhead;

% change the transformation matrix, such that it returns approximate SPM head coordinates
mri.transform     = itabvox2spmhead;
mri.vox2headOrig  = itabvox2itabhead;
mri.vox2head      = itabvox2spmhead;
mri.head2headOrig = spmhead2itabhead;
mri.coordsys      = 'spm';

% Do a second round of affine registration (rigid body) to get improved
% alignment with spm coordinate system. this is needed because there may be
% different conventions defining LPA and RPA. The affine registration may
% fail however, e.g. if the initial alignment is not close enough. In that
% case SPM will throw an error
if opt==1
  % use spm_affreg
  
  switch spm('ver')
      
    case 'SPM12'
      template = fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii');
        if ~exist('spm_affreg', 'file')
          addpath(fullfile(spm('Dir'),'toolbox','OldNorm'));
        end
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
  
  % the below, using just spm_affreg does not work robustly enough in some
  % cases
  flags.regtype = 'rigid';
  [M, scale]    = spm_affreg(V1,V2,flags);
  
  % some juggling around with the transformation matrices
  itabvox2spmhead2  = M \ V1.mat;
  spmhead2itabhead2 = itabvox2itabhead / itabvox2spmhead2;
   
  % update the transformation matrix
  mri.transform     = itabvox2spmhead2;
  
  % this one is unchanged
  mri.vox2headOrig  = itabvox2itabhead;
  
  % these are new
  mri.vox2head      = itabvox2spmhead2;
  mri.head2headOrig = spmhead2itabhead2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
  
elseif opt==2
  % use spm_normalise
  
  switch spm('ver')
    case 'SPM12'
      template = fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii');
        if ~exist('spm_affreg', 'file')
          addpath(fullfile(spm('Dir'),'toolbox','OldNorm'));
        end
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
  V1 = ft_write_mri(tname1, mri.anatomy,  'transform', mri.transform,  'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  V2 = ft_write_mri(tname2, mri2.anatomy, 'transform', mri2.transform, 'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  
  flags.nits       = 0; %set number of non-linear iterations to zero
  flags.regtype    = 'rigid';
  params           = spm_normalise(V2,V1,[],[],[],flags);
  spmhead2itabhead2 = spmhead2itabhead*V1.mat*params.Affine/V2.mat;
  itabvox2spmhead2  = spmhead2itabhead2\itabvox2itabhead;
  
  % update the transformation matrix
  mri.transform     = itabvox2spmhead2;
  
  % this one is unchanged
  mri.vox2headOrig  = itabvox2itabhead;
  
  % these are new
  mri.vox2head      = itabvox2spmhead2;
  mri.head2headOrig = spmhead2itabhead2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
end
