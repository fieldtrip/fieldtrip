function [warped]= sn2individual(P, input)

% SN2INDIVIDUAL warps the input coordinates (defined as Nx3 matrix) from
% normalised MNI coordinates to individual headspace coordinates, using the
% warp parameters defined in the structure spmparams.
%
% this is modified from code from nutmeg: nut_mni2mri, which was itself
% modified from code originally written by John Ashburner:
% http://www.sph.umich.edu/~nichols/JG2/get_orig_coord2.m

% Copyright (C) 2013-2017, Jan-Mathijs Schoffelen
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

  prm     = [3 3 3 0 0 0];
  Coef    = cell(1,3);
  Coef{1} = spm_bsplinc(P.Twarp(:,:,:,1),prm);
  Coef{2} = spm_bsplinc(P.Twarp(:,:,:,2),prm);
  Coef{3} = spm_bsplinc(P.Twarp(:,:,:,3),prm);
  
  VT        = P.tpm(1);
  M1        = VT.mat;
  d1        = VT.dim;
   
  M  = M1\P.Affine*P.image(1).mat; % this is an Affine mapping that goes from image voxels to TPM voxels
  d  = P.image(1).dim;
  
  [x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
  x3 = 1:d(3);
  
  y = zeros([d 3],'single');
  for z=1:length(x3)
    [t1,t2,t3] = defs(Coef,z,P.MT,prm,x1,x2,x3,M);
    tmp = M1(1,1)*t1 + M1(1,2)*t2 + M1(1,3)*t3 + M1(1,4);
    y(:,:,z,1) = tmp;
    tmp = M1(2,1)*t1 + M1(2,2)*t2 + M1(2,3)*t3 + M1(2,4);
    y(:,:,z,2) = tmp;
    tmp = M1(3,1)*t1 + M1(3,2)*t2 + M1(3,3)*t3 + M1(3,4);
    y(:,:,z,3) = tmp;
  end
  
  % 2b: write the deformation fields in x/y/z direction to temporary files
  V1.fname     = [tempname '.img'];
  V1.dim(1:3)  = P.image(1).dim(1:3);
  V1.pinfo     = [1 0 0]';
  V1.mat       = P.image(1).mat;
  V1.dt        = [64 0];
  V1.descrip   = 'Deformation field';
  spm_write_vol(V1,y(:,:,:,1));
  
  V2.fname     = [tempname '.img'];
  V2.dim(1:3)  = P.image(1).dim(1:3);
  V2.pinfo     = [1 0 0]';
  V2.mat       = P.image(1).mat;
  V2.dt        = [64 0];
  V2.descrip   = 'Deformation field';
  spm_write_vol(V2,y(:,:,:,2));
  
  V3.fname     = [tempname '.img'];
  V3.dim(1:3)  = P.image(1).dim(1:3);
  V3.pinfo     = [1 0 0]';
  V3.mat       = P.image(1).mat;
  V3.dt        = [64 0];
  V3.descrip   = 'Deformation field';
  spm_write_vol(V3,y(:,:,:,3));
 
  % first warp to voxel coordinates 
  warped = ft_warp_apply(inv(VT.mat), input);  % Express as voxel indices
  
  % apply the non-linear warp
  warped = cat(2, spm_sample_vol(V1,warped(:,1),warped(:,2),warped(:,3),1), ...
    spm_sample_vol(V2,warped(:,1),warped(:,2),warped(:,3),1), ...
    spm_sample_vol(V3,warped(:,1),warped(:,2),warped(:,3),1));
  
  % warp to head coordinates
  T      = P.image(1).mat;%*P.Affine;
  warped = ft_warp_apply(T, warped);
  
end

% the below is copied from spm_preproc_write8
%==========================================================================
% function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
%==========================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);

