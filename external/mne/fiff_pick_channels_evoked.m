function [res] = fiff_pick_channels_evoked(orig,include,exclude)
%
% [res] = fiff_pick_channels_evoked(orig,include,exclude)
%
% Pick desired channels from evoked-response data
%
% orig      - The original data
% include   - Channels to include (if empty, include all available)
% exclude   - Channels to exclude (if empty, do not exclude any)
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.9  2006/10/04 20:12:37  msh
%   Added fiff_read_evoked_all
%   Modified fiff_pick_channels_evoked to handle multiple data sets
%   Added bad channel writing to fiff_write_evoked
%
%   Revision 1.8  2006/06/29 22:12:28  msh
%   Fixed errors in channel picking
%
%   Revision 1.7  2006/04/26 00:14:24  msh
%   Return empty structure from fiff_pick_channels evoked and fiff_pick_types_evoked if there is no match.
%   Accuracy checking was incorrect in mne_add_coil_defs
%
%   Revision 1.6  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.5  2006/04/22 10:59:30  msh
%   Added fiff_pick_info
%
%   Revision 1.4  2006/04/20 23:33:20  msh
%   Separated fiff_pick_channels and fiff_pick_types from the _evoked versions
%
%   Revision 1.3  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.2  2006/04/14 15:49:49  msh
%   Improved the channel selection code and added ch_names to measurement info.
%
%   Revision 1.1  2006/04/14 00:45:42  msh
%   Added channel picking and fiff_invert_transform
%
%

me='MNE:fiff_pick_channels_evoked';

if nargin == 1
    res = orig;
    return;
elseif nargin == 2
    exclude = [];
elseif nargin ~= 3
    error(me,'Incorrect number of arguments');
end

if isempty(include) && isempty(exclude)
    res = orig;
    return;
end

sel = fiff_pick_channels(orig.info.ch_names,include,exclude);
if isempty(sel)
    res = [];
    fprintf(1,'Warning : No channels match the selection.\n');
    return;
end

res = orig;
%
%   Modify the measurement info
%
res.info = fiff_pick_info(res.info,sel);
%
%   Create the reduced data set
%
for k = 1:length(res.evoked)
    res.evoked(k).epochs = res.evoked(k).epochs(sel,:);
end

return;

end



