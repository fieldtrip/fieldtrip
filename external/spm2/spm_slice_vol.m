function varargout = spm_slice_vol(varargin)
% returns a slice through a memory mapped image - a compiled routine
% FORMAT X = spm_slice_vol(V,A,dim,hold);
% V      -  is a memory mapped image volume
% A      -  is a 4 x 4 transformation matrix
% dim    -  [i j] defines the two dimensions of the output image. The 
%           coordinates in 3-D space of the voxels in this image are assumed
%           to range from 1,1,0 to i,j,0.
% hold   -  sets the interpolation method for the resampling.
%           0          Zero-order hold (nearest neighbour).
%           1          First-order hold (trilinear interpolation).
%           2->127     Higher order Lagrange (polynomial) interpolation using
%                      different holds (second-order upwards).
%          -127 - -1   Different orders of sinc interpolation.
% X      -  output image
%_______________________________________________________________________
%
% spm_slice_vol returns a section through a memory mapped image volume
% on disk.  This section is the transverse slice at z = 0 after linear
% transformation according to matrix A
%
% see also spm_sample_vol.m
%
%_______________________________________________________________________
% @(#)spm_slice_vol.m	2.1 John Ashburner 99/04/19

%-This is merely the help file for the compiled routine
error('spm_slice_vol.c not compiled - see spm_MAKE.sh')
