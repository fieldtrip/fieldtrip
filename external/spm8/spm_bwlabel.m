function varargout = spm_bwlabel(varargin)
% Connected component labelling in 2D or 3D - a compiled routine
%
% FORMAT: [L,NUM] = spm_bwlabel(bw,conn)
% 
% The calling interface has been modelled on the 
% image processing toolbox routine bwlabel. See 
% the help for that routine for further information.
%
% Input:
% bw:         : (Binary) image to perform labelling on. Can
%               be 2 or 3D.
% conn        : Connectivity criterion. Could be 6(surface)
%               18(edge) or 26(corner). For a 2D bw these
%               correspond to 4, 8 and 8 respectively.
%
% Output:
% L           : Connected component image, i.e. image where
%               each non-zero voxel in bw will have a value
%               corresponding to its label.
% NUM         : Number of connected components in L.
%
%
% The implementation is not recursive (i.e. will no crash for
% large connected components) and is losely based on
% Thurfjell et al. 1992, A new three-dimensional connected
% components labeling algorithm with simultaneous object
% feature extraction capability. CVGIP: Graphical Models 
% and Image Processing 54(4):357-364.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson
% $Id$

%-This is merely the help file for the compiled routine
error('spm_bwlabel.c not compiled - see Makefile');




