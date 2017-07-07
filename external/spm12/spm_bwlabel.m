function [L,num] = spm_bwlabel(BW,n)
% Label connected components in 2D or 3D binary images - a compiled routine
%
% FORMAT [L,num] = spm_bwlabel(BW,n)
% BW     - 2D or 3D binary image to perform labelling on.
% n      - connectivity criterion: 6 (surface), 18 (edge) or 26 (corner).
%          [Default: 18].
%          (for a 2D image these correspond to 4, 8 and 8 respectively).
%
% L      - connected component image, i.e. image where each non-zero voxel
%          in BW will have a value corresponding to its label.
% num    - number of connected components in L.
%
%__________________________________________________________________________
%
% The implementation is loosely based on:
% Thurfjell et al. 1992, A new three-dimensional connected components
% labeling algorithm with simultaneous object feature extraction
% capability. CVGIP: Graphical Models and Image Processing 54(4):357-364.
%__________________________________________________________________________
% Copyright (C) 2002-2012 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id: spm_bwlabel.m 4929 2012-09-17 14:21:01Z guillaume $

%-This is merely the help file for the compiled routine
error('spm_bwlabel.c not compiled - see Makefile');
