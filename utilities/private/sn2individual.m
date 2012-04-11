function [warped]= sn2individual(P, input)

% SN2INDIVIDUAL warps the input coordinates (defined as Nx3 matrix) from
% normalised MNI coordinates to individual headspace coordinates, using the
% warp parameters defined in the structure spmparams.
%
% this is modified from code from nutmeg: nut_mni2mri, which was itself
% modified from code originally written by John Ashburner:
% http://www.sph.umich.edu/~nichols/JG2/get_orig_coord2.m

if numel(P.Tr)==0,
  % only an affine transformation has been done
  T      = P.VF.mat*P.Affine/(P.VG.mat);
  warped = warp_apply(T, input);
else
  % we need the spm_dctmtx function for the nonlinear case
  ft_hastoolbox('spm8', 1);

  dim  = P.VG.dim(1:3);
  xyz  = warp_apply(inv(P.VG.mat), input); % goes into voxel coordinates
  
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
  warped = warp_apply(T, xyz+xyztmp);
end;
