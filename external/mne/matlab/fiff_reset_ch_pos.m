function [res] = fiff_reset_ch_pos(chs)
%
% [res] = fiff_reset_ch_pos(chs)
%
% Reset channel position data to their original values as listed in
% the fif file
%
% NOTE: Only the coil_trans field is modified by this routine, not
% loc which remains to reflect the original data read from the fif file
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.4  2006/05/03 18:53:05  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.3  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.2  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.1  2006/04/13 22:37:03  msh
%   Added head_head_trans field to info.
%
%

me='MNE:fiff_reset_channel_pos';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin ~= 1
    error(me,'Wrong number of arguments');
end

res = chs;
for k = 1:length(res)
    loc = res(k).loc;
    if res(k).kind == FIFF.FIFFV_MEG_CH || ...
            res(k).kind == FIFF.FIFFV_REF_MEG_CH
        res(k).coil_trans = [ [ loc(4:6) loc(7:9) ...
            loc(10:12) ...
            loc(1:3) ] ; [ 0 0 0 1 ] ];
        res(k).coord_frame = FIFF.FIFFV_COORD_DEVICE;
    elseif res(k).kind == FIFF.FIFFV_EEG_CH
        res(k).eeg_loc = [ loc(1:3) loc(4:6) ];
        res(k).coord_frame = FIFF.FIFFV_COORD_HEAD;
    end
end

return;

end
