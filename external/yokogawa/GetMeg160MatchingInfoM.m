%  Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
% 
%  usage:
%    matching_info = GetMeg160MatchingInfoM( fid )
% 
%  arguments:
%    fid             : file ID
% 
%  return values:
%    matching_info        : Structure of MRI-MEG matching information.
%                        .meg_to_mri     : (4x4)transfer matrix of meg to mri
%                                            [xmri ymri zmri 1]' = meg_to_mri * [xmeg ymeg zmeg 1]', unit:[m]
%                        .mri_to_meg     : (4x4)transfer matrix of mri to meg 
%                                            [xmeg ymeg zmeg 1]' = meg_to_mri * [xmri ymri zmri 1]', unit:[m]
%                        .marker_count   : marker_count
%                        .marker         : 8 structures which has following infomation.
%                                            .mri_type, meg_type
%                                            .mri_done, meg_done
%                                            .mri_pos, meg_pos
%                        .marker_file_name
%  
%  confirmation of revision:
%   GetMeg160MatchingInfoM( Inf ) will show and return revision of this function.
%
