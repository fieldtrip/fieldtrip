function GX = spm_global(V)
% Compute the global mean for a volume image - a compiled routine
% FORMAT GX = spm_global(V)
% V   - image handle structure
% GX  - global mean
%__________________________________________________________________________
%
% spm_global returns the mean counts integrated over all the slices from
% the volume.
%
% The mean is estimated after discounting voxels outside the object using
% a criteria of greater than > (global mean)/8.
%__________________________________________________________________________
% Copyright (C) 1996-2012 Wellcome Trust Centre for Neuroimaging

% Anonymous
% $Id: spm_global.m 4921 2012-09-13 11:16:21Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_global.c not compiled - see Makefile')

% GX = zeros(numel(V),1);
% 
% for i=1:numel(V)
%     
%     D     = spm_data_read(V(i));
%     
%     iD    = isfinite(D);
%     
%     S     = mean(D(iD)) / 8;
%     
%     GX(i) = mean(D(iD & (D > S)));
%     
% end
