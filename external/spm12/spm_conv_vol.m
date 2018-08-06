function varargout = spm_conv_vol(varargin)
% Convolve a 3D volume with a three dimensional separable function
% FORMAT spm_conv_vol(V,Q,fx,fy,fz,offsets)
% V        -  input volume:
%             * a 3D array
%             * an image handle obtained by spm_vol
% Q        -  output volume:
%             * a 3D array (should probably be a lhs argument in this case)
%             * an image handle describing the format of the output image
% fx       -  the separable form of the function in x
% fy       -  the separable form of the function in y
% fz       -  the separable form of the function in z
% offsets  -  [i j k] contains the x, y and z shifts to reposition the
%             output
%__________________________________________________________________________
%
% spm_conv_vol is a compiled function (see spm_conv_vol.c).
%
% Separable means that f(x,y,z) = f(x)*f(y)*f(z) (= fx*fy*fz above)
%
% The convolution assumes zero padding in x and y with truncated smoothing 
% in z.
%
% If Q is an array with the same number of elements as the volume, the
% convolved volume will be stored there instead of on disk.  When Q 
% describes an output image, it is passed to the function spm_write_plane
% to write out each plane of the image.
%
% See also spm_conv.m and spm_smooth.m spm_write_plane.m
%
%__________________________________________________________________________
% Copyright (C) 1994-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner, Tom Nichols
% $Id: spm_conv_vol.m 6157 2014-09-05 18:17:54Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_conv_vol.c not compiled - see Makefile')
