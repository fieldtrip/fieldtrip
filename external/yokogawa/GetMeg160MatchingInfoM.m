% Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
%
% usage:
%   matching_info = GetMeg160MatchingInfoM( fid )
%
% arguments:
%   fid             : file ID
%
% return values:
%   matching_info        : Structure of MRI-MEG matching information.
%                       .meg_to_mri     : (4x4)transfer matrix of meg to mri
%                                           [xmri ymri zmri 1]' = meg_to_mri * [xmeg ymeg zmeg 1]', unit:[m]
%                       .mri_to_meg     : (4x4)transfer matrix of mri to meg 
%                                           [xmeg ymeg zmeg 1]' = meg_to_mri * [xmri ymri zmri 1]', unit:[m]
function    matching_info  =   GetMeg160MatchingInfoM( fid )