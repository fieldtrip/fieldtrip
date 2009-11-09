function GX = spm_global(V)
% returns the global mean for a memory mapped volume image - a compiled routine
% FORMAT GX = spm_global(V)
% V   - memory mapped volume
% GX  - mean global activity
%_______________________________________________________________________
%
% spm_global returns the mean counts integrated over all the  
% slices from the volume
%
% The mean is estimated after discounting voxels outside the object
% using a criteria of greater than > (global mean)/8
%
%_______________________________________________________________________
% @(#)spm_global.m	2.2 Anon 99/04/19

%-This is merely the help file for the compiled routine
error('spm_global.c not compiled - see spm_MAKE.sh')
