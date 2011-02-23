function [newdata] = mne_compensate_to(data,to)
%
% [newdata] = mne_compensate_to(data,to)
%
% Apply compensation to the data as desired
%


%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.4  2006/04/23 15:29:40  msh
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

me='MNE:mne_compensate_to';

newdata = data;
now     = mne_get_current_comp(newdata.info);
%
%   Are we there already?
%
if now == to
    fprintf(1,'Data are already compensated as desired\n');
end
%
%   Make the compensator and apply it to all data sets
%
comp = mne_make_compensator(newdata.info,now,to);
for k = 1:length(newdata.evoked)
   newdata.evoked(k).epochs = comp*newdata.evoked(k).epochs;
end
%
%  Update the compensation info in the channel descriptors
%
newdata.info.chs = mne_set_current_comp(newdata.info.chs,to);
return;

