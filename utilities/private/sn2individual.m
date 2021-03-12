function [warped]= sn2individual(P, input)

% SN2INDIVIDUAL warps the input coordinates (defined as Nx3 matrix) from
% normalised MNI coordinates to individual headspace coordinates, using the
% warp parameters defined in the structure spmparams.
%
% this is modified from code from nutmeg: nut_mni2mri, which was itself
% modified from code originally written by John Ashburner:
% http://www.sph.umich.edu/~nichols/JG2/get_orig_coord2.m

% Copyright (C) 2013-2021, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if isfield(P, 'Tr')
  % this is an old-style representation of the parameters, so it uses the
  % code adjusted from nut_mni2mri
  
  if numel(P.Tr)==0
    % only an affine transformation has been done
    T      = P.VF.mat*P.Affine/(P.VG.mat);
    warped = ft_warp_apply(T, input);
    
  else
    % we need the spm_dctmtx function for the nonlinear case
    if ~ft_hastoolbox('spm')
      % add SPM8 or later to the path
      ft_hastoolbox('spm8up', 1);
    end
    
    dim  = P.VG.dim(1:3);
    xyz  = ft_warp_apply(inv(P.VG.mat), input); % goes into voxel coordinates
    
    basX = spm_dctmtx(dim(1), size(P.Tr,1), xyz(:,1)-1);
    basY = spm_dctmtx(dim(2), size(P.Tr,2), xyz(:,2)-1);
    basZ = spm_dctmtx(dim(3), size(P.Tr,3), xyz(:,3)-1);
    
    siz = size(P.Tr);
    Tr1 = reshape(P.Tr(:,:,:,1),siz(1)*siz(2),siz(3));
    Tr2 = reshape(P.Tr(:,:,:,2),siz(1)*siz(2),siz(3));
    Tr3 = reshape(P.Tr(:,:,:,3),siz(1)*siz(2),siz(3));
    
    xyztmp = zeros(size(xyz));
    for i=1:size(xyz,1)
      bx = basX(i,:);
      by = basY(i,:);
      bz = basZ(i,:);
      tx = reshape(Tr1*bz', siz(1), siz(2) );
      ty = reshape(Tr2*bz', siz(1), siz(2) );
      tz = reshape(Tr3*bz', siz(1), siz(2) );
      xyztmp(i,:) = [bx*tx*by' bx*ty*by' bx*tz*by'];
    end
    
    T      = P.VF.mat*P.Affine;
    warped = ft_warp_apply(T, xyz+xyztmp);
  end
  
else
  % the only way I can come up with to do this, is to write a deformation
  % to disk, and to sample this one. This requires spm12 on the path
  ft_hastoolbox('spm12', 1);
   
  fprintf('creating the deformation field and writing it to a temporary file\n');
  
  fname  = [tempname,'.nii'];
  V      = nifti;
  V.dat  = file_array(fname, P.image.dim(1:3), [spm_type('float32') spm_platform('bigend')], 0, 1, 0);
  V.mat  = P.image.mat;
  if isfield(P.image, 'mat0') 
    V.mat0 = P.image.mat0;
  end
  V.descrip = 'deformation field';
  create(V);
  V.dat(:) = 0; % this is necessary, otherwise SPM fails: image too small
  
  P.image  = spm_vol(fname);
  spm_preproc_write8(P, zeros(6,4), [0 0], [0 1], 1, 1, nan(2,3), nan);
  
  [pth,nam,ext] = fileparts(fname);
  V  = nifti(fullfile(pth,['y_',nam,ext]));
  y  = squeeze(V.dat(:,:,:,:,:));
  
  siz = size(y);
  VT.dim = siz(1:3);
  VT.mat = P.tpm(1).mat;
  
  % 2b: write the deformation fields in x/y/z direction to temporary files
  V1.fname     = [tempname '.img'];
  V1.dim(1:3)  = VT.dim(1:3);
  V1.pinfo     = [1 0 0]';
  V1.mat       = VT.mat;
  V1.dt        = [64 0];
  V1.descrip   = 'Deformation field';
  spm_write_vol(V1,y(:,:,:,1));
  
  V2.fname     = [tempname '.img'];
  V2.dim(1:3)  = VT.dim(1:3);
  V2.pinfo     = [1 0 0]';
  V2.mat       = VT.mat;
  V2.dt        = [64 0];
  V2.descrip   = 'Deformation field';
  spm_write_vol(V2,y(:,:,:,2));
  
  V3.fname     = [tempname '.img'];
  V3.dim(1:3)  = VT.dim(1:3);
  V3.pinfo     = [1 0 0]';
  V3.mat       = VT.mat;
  V3.dt        = [64 0];
  V3.descrip   = 'Deformation field';
  spm_write_vol(V3,y(:,:,:,3));
 
  % first warp to voxel coordinates 
  input_vox = ft_warp_apply(inv(VT.mat), input);  % Express as voxel indices
  
  % apply the non-linear warp
  warped = cat(2, spm_sample_vol(V1,input_vox(:,1),input_vox(:,2),input_vox(:,3),1), ...
    spm_sample_vol(V2,input_vox(:,1),input_vox(:,2),input_vox(:,3),1), ...
    spm_sample_vol(V3,input_vox(:,1),input_vox(:,2),input_vox(:,3),1));
end
