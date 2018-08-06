function orig_coord = spm_get_orig_coord(coord, matname,PU)
% Determine corresponding co-ordinate in un-normalised image.
% FORMAT orig_coord = get_orig_coord2(coord, matname,PU)
% coord      - [x1 y1 z1 ; x2 y2 z2 ; etc] in MNI space (mm).
% matname    - File containing transformation information (_sn.mat).
%            - or the structure containing the transformation.
% PU         - Name of un-normalised image
% orig_coord - Co-ordinate in un-normalised image (voxel).
%
% FORMAT orig_coord = get_orig_coord2(coord, matname)
% coord      - [x1 y1 z1 ; x2 y2 z2 ; etc] in MNI space (mm).
% matname    - File containing transformation information (_sn.mat).
%            - or the structure containing the transformation.
% orig_coord - Original co-ordinate (mm).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_get_orig_coord.m 4873 2012-08-30 19:06:26Z john $

if ischar(matname)
    t = load(matname);
elseif isstruct(matname)
    t = matname;
else
    error('Wrong normalisation parameters!');
end
    
if size(coord,2)~=3, error('coord must be an N x 3 matrix'); end;
coord = coord';

Mat = inv(t.VG(1).mat);
xyz = Mat(1:3,:)*[coord ; ones(1,size(coord,2))];
Tr  = t.Tr;
Affine = t.Affine;
d   = t.VG(1).dim(1:3);

if nargin>2,
    VU   = spm_vol(PU);
    Mult = VU.mat\t.VF.mat*Affine;
%   disp('Output co-ordinates are in voxels');
else
    Mult = t.VF.mat*Affine;
%   disp('Output co-ordinates are in mm');
end;

if (prod(size(Tr)) == 0),
        affine_only = 1;
        basX = 0; tx = 0;
        basY = 0; ty = 0;
        basZ = 0; tz = 0;
else,
        affine_only = 0;
        basX = spm_dctmtx(d(1),size(Tr,1),xyz(1,:)-1);
        basY = spm_dctmtx(d(2),size(Tr,2),xyz(2,:)-1);
        basZ = spm_dctmtx(d(3),size(Tr,3),xyz(3,:)-1);
end;

if affine_only,
    xyz2 = Mult(1:3,:)*[xyz ; ones(1,size(xyz,2))];
else,
    for i=1:size(xyz,2),
        bx = basX(i,:);
        by = basY(i,:);
        bz = basZ(i,:);
        tx = reshape(...
            reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3))...
            *bz', size(Tr,1), size(Tr,2) );
        ty = reshape(...
            reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3))...
            *bz', size(Tr,1), size(Tr,2) );
        tz =  reshape(...
            reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3))...
            *bz', size(Tr,1), size(Tr,2) );
        xyz2(:,i) = Mult(1:3,:)*[xyz(:,i) + [bx*tx*by' ; bx*ty*by' ; bx*tz*by']; 1];
    end;
end;
orig_coord = xyz2';

return;
