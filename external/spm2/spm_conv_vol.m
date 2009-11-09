function varargout = spm_conv_vol(varargin)
% Convolves a mapped volume with a three dimensional separable function
% FORMAT spm_conv_vol(V,Q,fx,fy,fz,offsets)
% V        -  the input volume
%             - can be a 3D Matlab array, or an image handle
%               obtained by spm_vol
% Q        -  the output volume
%             - can be a 3D Matlab array (should probably be
%               a lhs argument in this case), or an image
%               handle describing the format of the output
%               image.
% fx       -  the separable form of the function in x
% fy       -  the separable form of the function in y
% fz       -  the separable form of the function in z
% offsets  - [i j k] contains the x, y and z shifts to reposition
%             the output
%_______________________________________________________________________
%
% spm_conv_vol is a compiled function (see spm_conv_vol.c).
%
% Separable means that f(x,y,z) = f(x)*f(y)*f(z) (= fx*fy*fz above)
%
% The convolution assumes zero padding in x and y with truncated smoothing 
% in z.
%
% If Q is an array with the same number of elements as the volume, the
% convolved volume will be store there instead of on disk.  When Q 
% describes an output image, it is passed to the function spm_write_plane
% to write out each plane of the image.
%
% See also spm_conv.m and spm_smooth.m spm_write_plane.m
%
%_______________________________________________________________________
% @(#)spm_conv_vol.m	2.1 John Ashburner, Tom Nichols 99/04/19

%-This is merely the help file for the compiled routine
error('spm_conv_vol.c not compiled - see spm_MAKE.sh')
