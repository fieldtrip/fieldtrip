function [Y1,Y2,Y3] = spm_invdef(X1,X2,X3,dimY,MX,MY)
% Create the inverse of a deformation field
% FORMAT [Y1,Y2,Y3] = spm_invdef(X1,X2,X3,dimY,MX,MY)
%   X1, X2 & X3 - Three components of original deformation field.
%                 Note that these point from voxels to a coordinate
%                 system in mm.
%   dimY        - A three element vector encoding the dimensions of
%                 the inverse of the deformation.
%   MX          - The voxel-to-world mapping for the elements of the
%                 original deformation.
%   MY          - The voxel-to-world mapping for the elements of the
%                 resulting inverse deformation.
%
%   Y1, Y2 & Y3 - Three components of inverse deformation field.
%                 Note that these point from voxels to a coordinate
%                 system in mm.
%
% Deformations are encoded as piecewise affine transforms. The space
% between each set of 8 neighbouring voxels is divided into 5
% tetrahedra, where there is an affine mapping within each of them.
%_______________________________________________________________________
% Inverting the deformation is as described in the appendix of:
%
%     John Ashburner, Jesper L.R. Andersson and Karl J. Friston.
%     "Image Registration Using a Symmetric Prior in Three Dimensions".
%     Human Brain Mapping 9:212-225(2000)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_invdef.m 4924 2012-09-13 16:55:09Z guillaume $

Y  = cat(4,X1,X2,X3);
Y  = spm_diffeo('invdef',Y,dimY,MX,MY);
Y1 = Y(:,:,:,1);
Y2 = Y(:,:,:,2);
Y3 = Y(:,:,:,3);

