function [mri] = align_neuromag2acpc(mri, method, template)

% ALIGN_NEUROMAG2ACPC performs an approximate alignment of the anatomical
% volume from NEUROMAG towards ACPC coordinates. Only the homogenous transformation
% matrix is modified and the coordsys-field is updated.
%
% Use as
%   mri = align_neuromag2acpc(mri)
%   mri = align_neuromag2acpc(mri, method)
%   mri = align_neuromag2acpc(mri, method, template)
%
% The first input argument is a FieldTrip MRI-structure, and the second optional
% argument specifies how the registration is to be done:
%   method = 0: only an approximate coregistration
%   method = 1: an approximate coregistration, followed by spm_affreg
%   method = 2: an approximate coregistration, followed by spm_normalise (default)
%
% When method = 1 or 2, an optional template filename can be specified, which denotes
% the filename of the target volume. This is required when running in deployed
% mode.
%
% See also ALIGN_CTF2ACPC, ALIGN_FSAVERAGE2MNI


if nargin<2
  method = 2;
end

%--------------------------------------------------------------------------
% do a first round of approximate coregistration using the predefined voxel
% locations of the fiducials and landmarks in the template T1 image from
% SPM

acpcvox2acpchead = [
  2     0     0   -92
  0     2     0  -128
  0     0     2   -74
  0     0     0     1
  ];

% these are the voxel indices of some points in the SPM canonical T1
acpcvox_Ac           = [46 64  37  1]';     % the anterior commissure
acpcvox_Ori          = [46 48  10  1]';     % approximately between the ears in T1.mnc
acpcvox_Nas          = [46 106 13  1]';     % approximately the nasion in T1.mnc
acpcvox_Lpa_canal    = [ 5 48  10  1]';     % Left ear canal
acpcvox_Rpa_canal    = [87 48  10  1]';     % Right ear canal

acpchead_Ac           = acpcvox2acpchead * acpcvox_Ac  ;
acpchead_Ori          = acpcvox2acpchead * acpcvox_Ori ;
acpchead_Nas          = acpcvox2acpchead * acpcvox_Nas ;
acpchead_Lpa_canal    = acpcvox2acpchead * acpcvox_Lpa_canal ;
acpchead_Rpa_canal    = acpcvox2acpchead * acpcvox_Rpa_canal ;

neuromagvox2neuromaghead  = mri.transform;
acpchead2neuromaghead     = ft_headcoordinates(acpchead_Nas(1:3), acpchead_Lpa_canal(1:3), acpchead_Rpa_canal(1:3), 'neuromag');

%neuromagvox2acpchead =  inv(acpchead2neuromaghead) * neuromagvox2neuromaghead;
neuromagvox2acpchead  =  acpchead2neuromaghead \ neuromagvox2neuromaghead;

% change the transformation matrix, such that it returns approximate SPM head coordinates
mri.transform     = neuromagvox2acpchead;
mri.vox2headOrig  = neuromagvox2neuromaghead;
mri.vox2head      = neuromagvox2acpchead;
mri.head2headOrig = acpchead2neuromaghead;
mri.coordsys      = 'acpc';

%--------------------------------------------------------------------------
% Do a second round of affine registration (rigid body) to get improved
% alignment with ACPC coordinate system. this is needed because there may be
% different conventions defining LPA and RPA. The affine registration may
% fail however, e.g. if the initial alignment is not close enough. In that
% case SPM will throw an error

if method==1
  % use spm_affreg
  
  switch lower(spm('ver'))
    case 'spm2'
      if isdeployed
        if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==1'); end
      else
        template = fullfile(spm('Dir'),'templates','T1.mnc');
      end
      
    case 'spm8'
      if isdeployed
        if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==1'); end
      else
        template = fullfile(spm('Dir'),'templates','T1.nii');
      end
      
    case 'spm12'
      if isdeployed
        if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==1'); end
      else
        template = fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii');
        if ~exist('spm_affreg', 'file')
          addpath(fullfile(spm('Dir'),'toolbox','OldNorm'));
        end
      end
      fprintf('using ''OldNorm'' affine registration\n');
      
    otherwise
      ft_error('unsupported SPM version');
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
  neuromagvox2acpchead2  = M \ V1.mat;
  acpchead2neuromaghead2 = neuromagvox2neuromaghead / neuromagvox2acpchead2;
  
  % update the transformation matrix
  mri.transform     = neuromagvox2acpchead2;
  
  % this one is unchanged
  mri.vox2headOrig  = neuromagvox2neuromaghead;
  
  % these are new
  mri.vox2head      = neuromagvox2acpchead2;
  mri.head2headOrig = acpchead2neuromaghead2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
  
elseif method==2
  % use spm_normalise
  
  switch lower(spm('ver'))
    case 'spm2'
      if isdeployed
        if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==2'); end
      else
        template = fullfile(spm('Dir'),'templates','T1.mnc');
      end
      
    case 'spm8'
      if isdeployed
        if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==2'); end
      else
        template = fullfile(spm('Dir'),'templates','T1.nii');
      end
      
    case 'spm12'
      % this uses the 'OldNorm' functionality, so the path needs to be
      % added, can only be done if non-deployed.
      if isdeployed
        if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==2'); end
      else
        template = fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii');
        if ~exist('spm_normalise', 'file')
          addpath(fullfile(spm('Dir'),'toolbox','OldNorm'));
        end
      end
      fprintf('using ''OldNorm'' normalisation\n');
      
    otherwise
      ft_error('unsupported SPM version');
  end
  mri2 = ft_read_mri(template);
  
  tname1 = [tempname, '.img'];
  tname2 = [tempname, '.img'];
  V1 = ft_write_mri(tname1, mri.anatomy,  'transform', mri.transform,  'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  V2 = ft_write_mri(tname2, mri2.anatomy, 'transform', mri2.transform, 'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  
  flags.nits       = 0; %set number of non-linear iterations to zero
  flags.regtype    = 'rigid';
  params           = spm_normalise(V2,V1,[],[],[],flags);
  acpchead2neuromaghead2 = acpchead2neuromaghead*V1.mat*params.Affine/V2.mat;
  neuromagvox2acpchead2  = acpchead2neuromaghead2\neuromagvox2neuromaghead;
  
  % update the transformation matrix
  mri.transform     = neuromagvox2acpchead2;
  
  % this one is unchanged
  mri.vox2headOrig  = neuromagvox2neuromaghead;
  
  % these are new
  mri.vox2head      = neuromagvox2acpchead2;
  mri.head2headOrig = acpchead2neuromaghead2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
end
