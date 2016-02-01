function [lf] = leadfield_fns(dip, vol)

% LEADFIELD_FNS calculates the FDM forward solution for a set of given
% dipolar sources
%
% [lf] = leadfield_fns(posin, vol, tol);
%
% with input arguments
%   dip    positions of the dipolar sources (MX3 matrix)
%   vol    structure of the volume conductor
%   tol    tolerance
% 
% The output argument lf 
% 
% The key elements of the vol structure are:
%   vol.condmatrix a 9XT (T tissues) matrix containing the conductivities
%   vol.seg        a segmented/labelled MRI
%   vol.deepelec   positions of the deep electrodes (NX3 matrix)
%    or
%   vol.bnd        positions of the external surface vertices

% The output lf is the forward solution matrix and contains the leadfields, calculated in 
% the N output voxels (vol.bnd or vol.deepelec) positions (a NXM matrix)

% Copyright (C) 2011, Cristiano Micheli and Hung Dang


% convert in voxel coordinates
dipvx = ft_warp_apply(inv(vol.transform),dip);
lf    = fns_leadfield(vol.transfer,vol.segdim,dipvx);
