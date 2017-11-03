function [warped]= individual2sn(P, input)

% INDIVIDUAL2SN warps the input coordinates (defined as Nx3 matrix) from
% individual headspace coordinates into normalised MNI coordinates, using the
% (inverse of the) warp parameters defined in the structure spmparams.
%
% this is code inspired by nutmeg and spm: nut_mri2mni, nut_spm_sn2def and
% nut_spm_invdef which were themselves modified from code originally written
% by John Ashburner:
% http://www.sph.umich.edu/~nichols/JohnsGems2.html
%
% Use as
%   [warped] = individual2sn(P, input)
%
% Input parameters:
%   P     = structure that contains the contents of an spm generated _sn.mat
%           file
%   input = Nx3 array containing the input positions

% Copyright (C) 2013, Jan-Mathijs Schoffelen
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
  % this is an old-style representation of the parameters, 
  
  % check for any version of SPM
  if ~ft_hastoolbox('spm')
    % add SPM8 to the path
    ft_hastoolbox('spm8', 1);
  end
  
  % The following is a three-step procedure
  
  % 1: create the spatial deformation field from the sn parameters
  [Def, M] = get_sn2def(P);
  
  % 2a: invert the spatial deformation field
  [p,f,e]         = fileparts(which('spm_invdef'));
  templatedirname = fullfile(p,'templates');
  d               = dir([templatedirname,'/T1*']);
  %VT              = spm_vol(fullfile(templatedirname,d(1).name));
  VT.dim = [91 109 91];
  VT.mat = [-2 0 0 92; 0 2 0 -128; 0 0 2 -74; 0 0 0 1];
  [iy(:,:,:,1), iy(:,:,:,2), iy(:,:,:,3)] = spm_invdef(Def{1},Def{2},Def{3},VT.dim(1:3),inv(VT.mat),M);
  
else
  
  % this requires spm12 on the path
  ft_hastoolbox('spm12', 1);
  
  prm     = [3 3 3 0 0 0];
  Coef    = cell(1,3);
  Coef{1} = spm_bsplinc(P.Twarp(:,:,:,1),prm);
  Coef{2} = spm_bsplinc(P.Twarp(:,:,:,2),prm);
  Coef{3} = spm_bsplinc(P.Twarp(:,:,:,3),prm);
  
  VT        = P.tpm(1);
  M1        = VT.mat;
  d1        = VT.dim;
   
  M  = M1\P.Affine*P.image(1).mat;
  d  = P.image(1).dim;
  
  [x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
  x3 = 1:d(3);
  
  y = zeros([d 3],'single');
  for z=1:length(x3)
    [t1,t2,t3] = defs(Coef,z,P.MT,prm,x1,x2,x3,M);
    y(:,:,z,1) = t1;
    y(:,:,z,2) = t2;
    y(:,:,z,3) = t3;
  end
  
  iy = spm_diffeo('invdef',y,d1,eye(4),P.image(1).mat);
  iy = spm_extrapolate_def(iy,M1);
  clear y;
end

% 2b: write the deformation fields in x/y/z direction to temporary files
V1.fname     = [tempname '.img'];
V1.dim(1:3)  = VT.dim(1:3);
V1.pinfo     = [1 0 0]';
V1.mat       = VT.mat;
V1.dt        = [64 0];
V1.descrip   = 'Inverse deformation field';
spm_write_vol(V1,iy(:,:,:,1));

V2.fname     = [tempname '.img'];
V2.dim(1:3)  = VT.dim(1:3);
V2.pinfo     = [1 0 0]';
V2.mat       = VT.mat;
V2.dt        = [64 0];
V2.descrip   = 'Inverse deformation field';
spm_write_vol(V2,iy(:,:,:,2));

V3.fname     = [tempname '.img'];
V3.dim(1:3)  = VT.dim(1:3);
V3.pinfo     = [1 0 0]';
V3.mat       = VT.mat;
V3.dt        = [64 0];
V3.descrip   = 'Inverse deformation field';
spm_write_vol(V3,iy(:,:,:,3));

% 3: extract the coordinates
warped = ft_warp_apply(inv(V1.mat), input);  % Express as voxel indices
warped = cat(2, spm_sample_vol(V1,warped(:,1),warped(:,2),warped(:,3),1), ...
  spm_sample_vol(V2,warped(:,1),warped(:,2),warped(:,3),1), ...
  spm_sample_vol(V3,warped(:,1),warped(:,2),warped(:,3),1));

%_______________________________________________________________________
function [Def,mat] = get_sn2def(sn)
% Convert a SPM _sn.mat file into a deformation field, and return it.
% This is code that was taken from SPM8.

dim = sn.VG(1).dim;
x   = 1:dim(1);
y   = 1:dim(2);
z   = 1:dim(3);
mat = sn.VG(1).mat;

[X,Y] = ndgrid(x,y);

st = size(sn.Tr);

if (prod(st) == 0)
    affine_only = true;
    basX = 0;
    basY = 0;
    basZ = 0;
else
    affine_only = false;
    basX = spm_dctmtx(sn.VG(1).dim(1),st(1),x-1);
    basY = spm_dctmtx(sn.VG(1).dim(2),st(2),y-1);
    basZ = spm_dctmtx(sn.VG(1).dim(3),st(3),z-1);
end,

Def = single(0);
Def(numel(x),numel(y),numel(z)) = 0;
Def = {Def; Def; Def};

for j=1:length(z)
    if (~affine_only)
        tx = reshape( reshape(sn.Tr(:,:,:,1),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
        ty = reshape( reshape(sn.Tr(:,:,:,2),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
        tz = reshape( reshape(sn.Tr(:,:,:,3),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );

        X1 = X    + basX*tx*basY';
        Y1 = Y    + basX*ty*basY';
        Z1 = z(j) + basX*tz*basY';
    end

    Mult = sn.VF.mat*sn.Affine;
    if (~affine_only)
        X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
        Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
        Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
    else
        X2= Mult(1,1)*X + Mult(1,2)*Y + (Mult(1,3)*z(j) + Mult(1,4));
        Y2= Mult(2,1)*X + Mult(2,2)*Y + (Mult(2,3)*z(j) + Mult(2,4));
        Z2= Mult(3,1)*X + Mult(3,2)*Y + (Mult(3,3)*z(j) + Mult(3,4));
    end

    Def{1}(:,:,j) = single(X2);
    Def{2}(:,:,j) = single(Y2);
    Def{3}(:,:,j) = single(Z2);
end;
%_______________________________________________________________________

% the below is copied from spm_preproc_write8, and adjusted a bit
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



              
              
