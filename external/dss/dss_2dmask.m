function [S] = dss_2dmask(X, mask, params)
% Perform masked DSS on 2d signals.
%  [S] = dss_2dmask(X, mask, params)
%    X       Mixed 2d signals as 3d array, 1st dimension as signal index
%    mask    2d denoising mask
%    params  Optional DSS algorithm parameters
%    S       Result 2d source signal estimates
%
%  Can be used eg. for spectrogram denoising.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

wdim = size(X,1)
fdim = size(X,2)
tdim = size(X,3)
X = reshape(X, wdim, fdim*tdim);
mask = reshape(mask, 1, fdim*tdim);

params.denf.h = @denoise_mask;
params.denf.params.mask = mask;

[state, A, W] = denss(X, params);

S = W * X;

sdim = size(S,1);

S = reshape(sdim, fdim, tdim);
