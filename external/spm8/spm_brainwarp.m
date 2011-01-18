function varargout = spm_brainwarp(varargin)
% Part of nonlinear spatial normalisation - a compiled routine
%_______________________________________________________________________
% [Alpha,Beta,Var] = spm_brainwarp(VG,VF,Affine,basX,basY,basZ,...
%                                   dbasX,dbasY,dbasZ,T,fwhm,VW,VW2)
% VG    - Template volume(s) (see spm_vol)
% VF    - Object volume
% Affine    - The affine transformation which maps between the object
%     and template.
% basX  - Basis vectors in X. # rows must eq. VG(1)
% basY  - Basis vectors in Y. # rows must eq. VG(2)
% basZ  - Basis vectors in Z. # rows must eq. VG(3)
% dbasX - Derivatives of basis vectors in X. # rows must eq. VG(1)
% dbasY - Derivatives of basis vectors in Y. # rows must eq. VG(2)
% dbasZ - Derivatives of basis vectors in Z. # rows must eq. VG(3)
% T - The current parameter estimates.
% fwhm  - The approximate smoothness of the images.
% VW    - an optional weighting volume for determining which voxels
%         should be weighted more heavily in the fitting process.
%         This volume should have the same dimensions and position
%         as the volumes in VG.
% VW2   - another optional weighting volume for determining which voxels
%         should be weighted more heavily in the fitting process.
%         This volume should have the same dimensions and position
%         as the volumes in VF.
% Without the weighting volumes, all voxels are assigned weights that
% are uniformly one.
% 
% Alpha - A*A - where A is the design matrix
% Beta  - A*b - where f is the object image
% Var   - the approximate chi^2 (corrected for number of resels).
%_______________________________________________________________________
% 
% The voxels of g1, g2.. are sampled according to the smoothness of the
% image (fwhm). The corresponding voxels of f are determined according
% to the current parameter estimates and the affine transform.  See
% "spm_write_sn.m" for more details about how this is done.
% 
% 
%-----------------------------------------------------------------------
% 
% The design matrix A is generated internally from:
% 
% diag(w)*[diag(df/dx)*B diag(df/dy)*B diag(df/dz)*B ...
%      diag(g1)*[1 x y z] ...
%      diag(g2)*[1 x y z] ...]
% 
% where df/dx, df/dy & df/dz are column vectors containing the gradient
%   of image f with respect to displacements in x, y & z
%   (in the space of g).
% 
%   B is generated from kron(basZ,kron(basY,BasX)). Each column of
%   B is a basis image.
% 
%   g1, g2.. are template images.
% 
%   x, y & z are simply the spatial coordinates of the voxels of f.
% 
%   s1, s2.. are the current estimates for the required scaling
%   factors. These are derived from T(3*prod(VG(1:3))+1),
%   T(3*prod(VG(1:3))+2)...
%
%       w is an optional vector of weights, where w = 1/(1/w1 + 1/w2)
%       where w1 and w2 are derived from the optional weighting images.
% 
% The vector b contains [diag(w)*(f - diag(g1)*s1 - diag(g1)*x*s2 - ...)].
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_brainwarp.m 1143 2008-02-07 19:33:33Z spm $


%-This is merely the help file for the compiled routine
error('spm_brainwarp.c not compiled - see Makefile')
