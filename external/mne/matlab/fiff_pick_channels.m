function [sel] = fiff_pick_channels(ch_names,include,exclude)
%
% [sel] = fiff_pick_channels(ch_names,include,exclude)
%
% Make a selector to pick desired channels from data
%
% ch_names  - The channel name list to consult
% include   - Channels to include (if empty, include all available)
% exclude   - Channels to exclude (if empty, do not exclude any)
%
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.5  2008/04/06 16:35:35  msh
%   Added reasonable handling of duplicate channel names.
%
%   Revision 1.4  2006/06/29 22:12:28  msh
%   Fixed errors in channel picking
%
%   Revision 1.3  2006/04/27 22:38:37  msh
%   Splitting an empty list now results in an empty output.
%   Added fiff_write_double and fiff_write_short
%   Write an array of floats, ints, and doubles instead of just one value.
%   Fixed help text of fiff_pick_channels.
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/20 23:33:20  msh
%   Separated fiff_pick_channels and fiff_pick_types from the _evoked versions
%
%

me='MNE:fiff_pick_channels';

nchan = length(ch_names);
if nargin == 1
    sel = ones(1,nchan);
    for k = 1:nchan
        sel(k) = k;
    end
    return;
elseif nargin == 2
    exclude = [];
elseif nargin ~= 3
    error(me,'Incorrect number of arguments');
end

if isempty(include)
    %
    %   Include all initially
    %
    sel = zeros(1,nchan);
    for k = 1:nchan
        sel(k) = k;
    end
    nzero = 0;
    for k = 1:length(exclude)
        c = strmatch(exclude{k},ch_names,'exact');
        nzero = 0;
        if length(c) > 0
            sel(c(1)) = 0;
            nzero = nzero + 1;
        end
    end
    %
    %  Check for exclusions
    %
    if nzero > 0
        newsel = zeros(1,nchan-nzero);
        p = 0;
        for k = 1:nchan
            if sel(k) > 0
                p = p + 1;
                newsel(p) = sel(k);
            end
        end
        sel = newsel;
    end
else
    %
    %   First do the channels to be included
    %
    sel = zeros(1,length(include));
    nzero = 0;
    for k = 1:length(include)
        c = strmatch(include{k},ch_names,'exact');
        if ~length(c)
            error(me,'Missing channel %s',include{k});
        elseif length(c) > 1
            disp(sprintf('Ambiguous channel, taking first occurence: %s',include{k}));
        end
        %
        %  Is this channel in the exclusion list?
        %
        sel(k) = c(1);
        if ~isempty(exclude)
            c = strmatch(include{k},exclude,'exact');
            if length(c) > 0
                sel(k) = 0;
                nzero = nzero + 1;
            end
        end
    end
    %
    %    Check whether some channels were excluded
    %
    if nzero > 0
        newsel = zeros(1,length(include)-nzero);
        p = 0;
        for k = 1:length(include)
            if sel(k) > 0
                p = p + 1;
                newsel(p) = sel(k);
            end
        end
        sel = newsel;
    end
end

return;

end



