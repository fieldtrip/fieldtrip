function [hemi] = mne_find_source_space_hemi(src)
%
% function mne_find_source_space_hemi(src)
%
% Return the hemisphere id for a source space
%
% src      - The source space to investigate
% hemi     - Deduced hemisphere id
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.2  2006/09/25 19:48:16  msh
%   Added projection item kinds to fiff_define_constants
%   Changed some fields to int32 in surface structures
%
%   Revision 1.1  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%

me='MNE:mne_find_source_space_hemi';

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

if nargin ~= 1
    error(me,'Incorrect number of arguments');
end

xave = sum(src.rr(:,1));

if xave < 0
    hemi = int32(FIFF.FIFFV_MNE_SURF_LEFT_HEMI);
else
    hemi = int32(FIFF.FIFFV_MNE_SURF_RIGHT_HEMI);
end

return;

