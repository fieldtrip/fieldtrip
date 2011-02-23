function [sel] = fiff_pick_types(info,meg,eeg,stim,include,exclude)
%
% [sel] = fiff_pick_types(info,meg,eeg,stim,exclude)
%
% Create a selector to pick desired channel types from data
%
% info      - The measurement info
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
%   Revision 1.5  2006/06/29 22:12:28  msh
%   Fixed errors in channel picking
%
%   Revision 1.4  2006/05/03 18:53:04  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.3  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.2  2006/04/21 11:33:17  msh
%   Added raw segment reading routines
%
%   Revision 1.1  2006/04/20 23:33:20  msh
%   Separated fiff_pick_channels and fiff_pick_types from the _evoked versions
%
%

me='MNE:fiff_pick_types';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

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

pick = zeros(1,info.nchan);

for k = 1:info.nchan
    kind = info.chs(k).kind;
    if (kind == FIFF.FIFFV_MEG_CH || kind == FIFF.FIFFV_REF_MEG_CH) && meg
        pick(k) = true;
    elseif kind == FIFF.FIFFV_EEG_CH && eeg
        pick(k) = true;
    elseif kind == FIFF.FIFFV_STIM_CH && stim
        pick(k) = true;
    end
end

p = 0;
for k = 1:info.nchan
    if pick(k)
        p = p + 1;
        myinclude{p} = info.ch_names{k};
    end
end

if ~isempty(include)
    for k = 1:length(include)
        p = p + 1;
        myinclude{p} = include{k};
    end
end

if p == 0
    sel = [];
else
    sel = fiff_pick_channels(info.ch_names,myinclude,exclude);
end

return;

end
