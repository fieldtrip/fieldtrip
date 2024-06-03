function flip = spm_flip_analyze_images
% Do Analyze format images need to be left-right flipped? The default
% behaviour is to have the indices of the voxels stored as left-handed and
% interpret the mm coordinates within a right-handed coordinate system.
%
% Note: the behaviour used to be set in spm_defaults.m, but this has now
% been changed.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2002-2022 Wellcome Centre for Human Neuroimaging

flip = 1;
