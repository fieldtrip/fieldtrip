function varargout = spm_sample_vol(varargin)
% returns voxel values from a memory mapped image - a compiled routine
% FORMAT X = spm_sample_vol(V,x,y,z,hold);
% V      -  is a memory mapped image volume
% x      -  matrix of x coordinates {pixels}
% y      -  matrix of y coordinates {pixels}
% z      -  matrix of z coordinates {pixels}
% hold   -  sets the interpolation method for the resampling.
%           0          Zero-order hold (nearest neighbour).
%           1          First-order hold (trilinear interpolation).
%           2->127     Higher order Lagrange (polynomial) interpolation using
%                      different holds (second-order upwards). 
%          -127 - -1   Different orders of sinc interpolation. 
% X      -  output image
%
% OR     [X,dX,dY,dZ] = spm_sample_vol(V,x,y,z,hold);
% Similar to above, except that the derivatives in the three orthogonal
% directions are also returned.
%_______________________________________________________________________
%
% spm_sample_vol will return the voxel values from a memory mapped
% volume indicated by V at coordinates x,y,z.  Values from coordinates
% outside the image are set to zero. x, y and z must be matrices of the
% same dimensions
%
% see also spm_slice_vol.m
%
%_______________________________________________________________________
% @(#)spm_sample_vol.m	2.1 John Ashburner 99/04/19

%-This is merely the help file for the compiled routine
error('spm_sample_vol.c not compiled - see spm_MAKE.sh')
