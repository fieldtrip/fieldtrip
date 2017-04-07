function varargout = spm_slice_vol(varargin)
% Return a section through an image volume - a compiled routine
% FORMAT X = spm_slice_vol(V,A,dim,hold)
% V        -  spm_vol structure
% A        -  4 x 4 transformation matrix
% dim      -  [i j] defining the two dimensions of the output image.
%             The coordinates in 3-D space of the voxels in this image are
%             assumed to range from 1,1,0 to i,j,0
% hold     -  interpolation method for the resampling:
%              0         : Zero-order hold (nearest neighbour)
%              1         : First-order hold (trilinear interpolation)
%              2->127    : Higher order Lagrange (polynomial) interpolation
%                          using different holds (second-order upwards)
%              -127 - -1 : Different orders of sinc interpolation
%
% X        -  output image
%__________________________________________________________________________
%
% spm_slice_vol returns a section through an image volume.
% This section is the transverse slice at z = 0 after linear transformation
% according to matrix A
%
% See also: spm_vol.m, spm_sample_vol.m
%__________________________________________________________________________
% Copyright (C) 1999-2015 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_slice_vol.m 6340 2015-02-16 12:25:56Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_slice_vol.c not compiled - see Makefile')
