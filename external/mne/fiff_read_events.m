function [events, mappings] = fiff_read_events(source, tree)
%
% [events, mappings] = fiff_read_events(source, tree)
%
% Read the events
%
% If tree is specified, source is assumed to be an open file id,
% otherwise a the name of the file to read. If tree is missing, the
% meas output argument should not be specified.
%
%
%   Author : Martin Luessi, MGH Martinos Center
%   License : BSD 3-clause


global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_read_events';

if nargin ~= 2 & nargin ~= 1
    error(me,'Incorrect number of arguments');
end

if nargin == 1
    [ fid, tree ] = fiff_open(source);
    open_here = true;
else
    fid = source;
    open_here = false;
end

%
%   Find the desired blocks
%
data = fiff_dir_tree_find(tree, FIFF.FIFFB_MNE_EVENTS);
if length(data) == 0
    if open_here
        fclose(fid);
    end
    error(me,'Could not find events');
end

%
%   Read the events
%
for k = 1:data.nent
    kind = data.dir(k).kind;
    pos  = data.dir(k).pos;
    if kind == FIFF.FIFF_MNE_EVENT_LIST
        tag = fiff_read_tag(fid, pos);
        events = tag.data;
        break
    end
end

if ~exist('events','var')
    if open_here
        fclose(fid);
    end
    error(me,'Events not found');
end

events = reshape(events, 3, length(events) / 3)';

%
% Read the event mappings
%
mappings = '';
data = fiff_dir_tree_find(tree, FIFF.FIFFB_MNE_EVENTS);

for k = 1:data.nent
    kind = data.dir(k).kind;
    pos  = data.dir(k).pos;
    if kind == FIFF.FIFF_DESCRIPTION
        tag = fiff_read_tag(fid, pos);
        mappings = tag.data;
        break
    end
end

if open_here
    fclose(fid);
end

return;



