function [mri] = align_ctf2spm(mri, opt, template)

% ALIGN_CTF2SPM performs an approximate alignment of the anatomical volume from CTF
% towards SPM coordinates. Only the homogeneous transformation matrix is modified and
% the coordsys-field is updated.
%
% Use as
%   mri = align_ctf2spm(mri)
%   mri = align_ctf2spm(mri, opt)
%   mri = align_ctf2spm(mri, opt, template)
% where mri is a FieldTrip MRI-structure, and opt is an optional argument specifying
% how the registration is to be done:
%   opt = 0: only an approximate coregistration
%   opt = 1: an approximate coregistration, followed by spm_affreg
%   opt = 2 (default): an approximate coregistration, followed by spm_normalise
%
% When opt = 1 or 2, an optional template filename can be specified, which denotes
% the filename of the target volume. This option is required when running in deployed
% mode.

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

ctfvox2ctfhead  = mri.transform;
spmhead2ctfhead = ft_headcoordinates(spmhead_Nas(1:3), spmhead_Lpa_canal(1:3), spmhead_Rpa_canal(1:3), 'ctf');

%ctfvox2spmhead =  inv(spmhead2ctfhead) *  ctfvox2ctfhead;
ctfvox2spmhead =  spmhead2ctfhead \ ctfvox2ctfhead;

% change the transformation matrix, such that it returns approximate SPM head coordinates
mri.transform     = ctfvox2spmhead;
mri.vox2headOrig  = ctfvox2ctfhead;
mri.vox2head      = ctfvox2spmhead;
mri.head2headOrig = spmhead2ctfhead;
mri.coordsys      = 'spm';

% Do a second round of affine registration (rigid body) to get improved
% alignment with spm coordinate system. this is needed because there may be
% different conventions defining LPA and RPA. The affine registration may
% fail however, e.g. if the initial alignment is not close enough. In that
% case SPM will throw an error

% check for any version of SPM
if ~ft_hastoolbox('spm')
  % add SPM8 to the path
  ft_hastoolbox('spm8', 1);
end


if opt==1
  % use spm_affreg
  
  switch lower(spm('ver'))
    case 'spm2'
      if isdeployed
        if nargin<3, error('you need to specify a template filename when in deployed mode and using opt==1'); end
      else
        template = fullfile(spm('Dir'),'templates','T1.mnc');
      end
      
    case 'spm8'
      if isdeployed
        if nargin<3, error('you need to specify a template filename when in deployed mode and using opt==1'); end
      else
        template = fullfile(spm('Dir'),'templates','T1.nii');
      end
      
    case 'spm12'
      if isdeployed
        if nargin<3, error('you need to specify a template filename when in deployed mode and using opt==1'); end
      else
        template = fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii');
        if ~exist('spm_affreg', 'file')
          addpath(fullfile(spm('Dir'),'toolbox','OldNorm'));
        end
      end
      fprintf('using ''OldNorm'' affine registration\n');

    otherwise
      error('unsupported spm-version');
  end
  mri2 = ft_read_mri(template);
  
  tname1 = [tempname, '.img'];
  tname2 = [tempname, '.img'];
  V1 = ft_write_mri(tname1, mri.anatomy,  'transform', mri.transform,  'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  V2 = ft_write_mri(tname2, mri2.anatomy, 'transform', mri2.transform, 'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  
  % the below, using just spm_affreg does not work robustly enough in some cases
  flags.regtype = 'rigid';
  [M, scale]    = spm_affreg(V1,V2,flags);
  
  % some juggling around with the transformation matrices
  ctfvox2spmhead2  = M \ V1.mat;
  spmhead2ctfhead2 = ctfvox2ctfhead / ctfvox2spmhead2;
  
  % update the transformation matrix
  mri.transform     = ctfvox2spmhead2;
  
  % this one is unchanged
  mri.vox2headOrig  = ctfvox2ctfhead;
  
  % these are new
  mri.vox2head      = ctfvox2spmhead2;
  mri.head2headOrig = spmhead2ctfhead2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
  
elseif opt==2
  % use spm_normalise
  
  switch lower(spm('ver'))
    case 'spm2'
      if isdeployed
        if nargin<3, error('you need to specify a template filename when in deployed mode and using opt==2'); end
      else
        template = fullfile(spm('Dir'),'templates','T1.mnc');
      end
      
    case 'spm8'
      if isdeployed
        if nargin<3, error('you need to specify a template filename when in deployed mode and using opt==2'); end
      else
        template = fullfile(spm('Dir'),'templates','T1.nii');
      end
      
    case 'spm12'
      % this uses the 'OldNorm' functionality, so the path needs to be
      % added, can only be done if non-deployed.
      if isdeployed
        if nargin<3, error('you need to specify a template filename when in deployed mode and using opt==2'); end
      else
        template = fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii');
        if ~exist('spm_normalise', 'file')
          addpath(fullfile(spm('Dir'),'toolbox','OldNorm'));
        end
      end
      fprintf('using ''OldNorm'' normalisation\n');
      
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
  spmhead2ctfhead2 = spmhead2ctfhead*V1.mat*params.Affine/V2.mat;
  ctfvox2spmhead2  = spmhead2ctfhead2\ctfvox2ctfhead;
  
  % update the transformation matrix
  mri.transform     = ctfvox2spmhead2;
  
  % this one is unchanged
  mri.vox2headOrig  = ctfvox2ctfhead;
  
  % these are new
  mri.vox2head      = ctfvox2spmhead2;
  mri.head2headOrig = spmhead2ctfhead2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
end
