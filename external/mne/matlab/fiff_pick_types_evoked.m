function [res] = fiff_pick_types_evoked(orig,meg,eeg,stim,include,exclude)
%
% [res] = fiff_pick_types_evoked(orig,meg,eeg,stim,include,exclude)
%
% Pick desired channels types from evoked-response data
%
% orig      - The original data
% meg       - Include MEG channels
% eeg       - Include EEG channels
% stim      - Include stimulus channels
% include   - Additional channels to include (if empty, do not add any)
% exclude   - Channels to exclude (if empty, do not exclude any)
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.10  2008/06/12 19:59:19  msh
%   Fixed the help text
%
%   Revision 1.9  2007/12/23 17:01:24  msh
%   Fixed error handling multiple datasets
%
%   Revision 1.8  2006/04/26 00:14:24  msh
%   Return empty structure from fiff_pick_channels evoked and fiff_pick_types_evoked if there is no match.
%   Accuracy checking was incorrect in mne_add_coil_defs
%
%   Revision 1.7  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.6  2006/04/21 11:33:17  msh
%   Added raw segment reading routines
%
%   Revision 1.5  2006/04/20 23:33:20  msh
%   Separated fiff_pick_channels and fiff_pick_types from the _evoked versions
%
%   Revision 1.4  2006/04/17 11:52:15  msh
%   Added coil definition stuff
%
%   Revision 1.3  2006/04/15 12:21:00  msh
%   Several small improvements
%
%   Revision 1.2  2006/04/14 15:49:49  msh
%   Improved the channel selection code and added ch_names to measurement info.
%
%   Revision 1.1  2006/04/14 00:45:42  msh
%   Added channel picking and fiff_invert_transform
%
%

me='MNE:fiff_pick_types_evoked';

if nargin == 5
    exclude = [];
elseif nargin == 4
    include = [];
    exclude = [];
elseif nargin == 3
    include = [];
    exclude = [];
    stim    = false;
elseif nargin == 2
    include = [];
    exclude = [];
    stim    = false;
    eeg     = false;
elseif nargin ~= 6
    error(me,'Incorrect number of arguments');
end

sel = fiff_pick_types(orig.info,meg,eeg,stim,include,exclude);

if isempty(sel)
    res = [];
    fprintf(1,'Warning : No channels match the selection.\n');
    return;
end

res = orig;
res.info.chs      = res.info.chs(sel);
res.info.ch_names = res.info.ch_names(sel);
%
%   Create the reduced data set
%
for k = 1:length(res.evoked)
    res.evoked(k).epochs = res.evoked(k).epochs(sel,:);
end
res.info.nchan    = length(sel);

return;

end
