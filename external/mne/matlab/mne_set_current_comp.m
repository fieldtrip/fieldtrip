function new_chs = mne_set_current_comp(chs,value)
%
% mne_set_current_comp(chs,value)
%
% Set the current compensation value in the channel info structures
%


%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.5  2008/11/18 02:38:51  msh
%   Modified mne_ex_read_evoked to apply projection and compensation
%   Modified mne_ex_read_raw to call mne_set_current_comp
%
%   Revision 1.4  2006/04/23 15:29:41  msh
%   Added MGH to the copyright
%
%   Revision 1.3  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.2  2006/04/14 15:49:49  msh
%   Improved the channel selection code and added ch_names to measurement info.
%
%   Revision 1.1  2006/04/12 10:51:19  msh
%   Added projection writing and compensation routines
%
%

me='MNE:mne_set_current_comp';

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end

new_chs = chs;

lower_half = hex2dec('FFFF');
for k = 1:length(chs)
    if chs(k).kind == FIFF.FIFFV_MEG_CH
        coil_type = bitand(double(chs(k).coil_type),lower_half);
        new_chs(k).coil_type = int32(bitor(coil_type,bitshift(value,16)));
    end
end

return;

