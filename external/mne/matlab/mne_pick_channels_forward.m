function [fwd] = mne_pick_channels_forward(orig,include,exclude)
%
% [fwd] = mne_pick_channels_forward(orig,include,exclude)
%
% Pick desired channels from a forward solution
%
% orig      - The original forward solution
% include   - Channels to include (if empty, include all available)
% exclude   - Channels to exclude (if empty, do not exclude any)
%
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.2  2009/03/04 22:04:31  msh
%   Fixed typos in comments
%
%   Revision 1.1  2009/02/11 22:04:02  msh
%   New routine to pick selected channels from the forward solution
%
%

me='MNE:mne_pick_channels_forward';

if nargin == 1
    fwd = orig;
    return;
elseif nargin == 2
    exclude = [];
elseif nargin ~= 3
    error(me,'Incorrect number of arguments');
end

if isempty(include) && isempty(exclude)
    fwd = orig;
    return;
end

if isempty(orig.sol.row_names) 
    error(me,'Cannot pick from a forward solution without channel names');
end

fwd = orig;
%
%   First do the channels to be included
%
sel = fiff_pick_channels(orig.sol.row_names, include, exclude);
if isempty(sel)
   error(me, 'Nothing remains after picking');
end
fwd.sol.data = fwd.sol.data(sel,:);
fwd.chs = fwd.chs(sel);
%
%   Select the desired stuff
%
for p = 1:length(sel)
   names{p} = fwd.sol.row_names{sel(p)};
end
fwd.sol.row_names = names;
fwd.sol.nrow      = length(sel);
fwd.nchan         = length(sel);
%
return;
