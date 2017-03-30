function varargout = spm_mrf(varargin)
% Markov Random Field Code - a compiled routine
%_______________________________________________________________________
%
% FORMAT q1 = spm_mrf(q0,p,G,w)
% q0 - Original responsibilities.
%      Values are stored as uint8 and converted to responsibilities by
%      rescaling by 1/255. Dimensions are dim1 x dim2 x dim3 x K
% p  - Probabilities.
%      Values are as single precision floating point.  This array must
%      have the same dimensions as q0.
% G  - Matrix used to encode neighbourhood information.
%      May be of a number of types.
%        i) K x K matrix (single precision), where K is the 4th
%           dimension of q0 and p.  This matrix is shared by all voxels.
%       ii) K x 1 vector (single precision), encoding the diagonal of
%           a matrix.
%      iii) dim1 x dim2 x dim3 x K x K (single precision).  Encodes a
%           different matrix at each voxel, and is very memory hungry.
%       iv) dim1 x dim2 x dim3 x (K*(K-1)/2) (single precision).
%           Encodes a symmetric matrix, with zeros on the diagonal, at
%           each voxel.
%        v) dim1 x dim2 x dim3 x (K*(K-1)/2) (uint8).
%           Encodes a symmetric matrix, with zeros on the diagonal, at
%           each voxel. Saves more memory by using uint8.  Note that
%           when used, the uint8 values are rescaled by -1/(2^4).
% w  - A vector of three weights, which normally encode the reciprocal
%      of the square of the voxel sizes.  This is for dealing with
%      anisotropic voxels. If this argument is not supplied, then
%      [1 1 1] is assumed.
% q1 - Output responsibilities.
%
% FORMAT spm_mrf(q,p,G,w)
% This is the dodgy way of using the function, as it changes the RHS
% argument (q) and can lead to some strange side effects.  This approach
% would not be ebdorsed by the MathWorks, but it does save a bit of memory.
%
%
% The MRF updates are done using a red-black checkerboard scheme.  Each
% voxel is updated by q = (exp(G'*a).*p)/sum(exp(G'*a).*p), where
% vector a is computed from the number of neighbours of each type
% (divided by 6). The contribution of each neighbour is scaled by w.
%
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_mrf.m 6881 2016-09-19 09:48:54Z john $

%-This is merely the help file for the compiled routine
error('spm_mrf.c not compiled - see Makefile')
