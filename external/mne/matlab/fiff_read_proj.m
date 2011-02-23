function [ projdata ] = fiff_read_proj(fid,node)

%
% [ projdata ] = fiff_read_proj(fid,node)
%
% Read the SSP data under a given directory node
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.7  2006/05/16 00:39:32  msh
%   Fixed error in mne_ex_read_raw: All projection items were not activated.
%   List the initial states of projection items when they are loaded.
%
%   Revision 1.6  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.5  2006/04/21 15:06:03  msh
%   Do not list projection item status because it might be confusing.
%
%   Revision 1.4  2006/04/21 14:43:59  msh
%   Report projection operators.
%   Some more checks in raw data reading.
%
%   Revision 1.3  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.2  2006/04/12 10:29:02  msh
%   Made evoked data writing compatible with the structures returned in reading.
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_read_proj';

if nargin ~= 2
    error(me,'Incorrect number of arguments');
end

projdata = struct('kind',{},'active',{},'desc',{},'data',{});
%
%   Locate the projection data
%
nodes = fiff_dir_tree_find(node,FIFF.FIFFB_PROJ);
if length(nodes) == 0
    return;
end
tag = find_tag(nodes(1),FIFF.FIFF_NCHAN);
if ~isempty(tag)
    global_nchan = tag.data;
end
items = fiff_dir_tree_find(nodes(1),FIFF.FIFFB_PROJ_ITEM);
for i = 1:length(items)
    %
    %   Find all desired tags in one item
    %
    item = items(i);
    tag = find_tag(item,FIFF.FIFF_NCHAN);
    if ~isempty(tag)
        nchan = tag.data;
    else
        nchan = global_nchan;
    end
    tag = find_tag(item,FIFF.FIFF_DESCRIPTION);
    if ~isempty(tag)
        desc = tag.data;
    else
        tag = find_tag(item,FIFF.FIFF_NAME);
        if ~isempty(tag)
            desc = tag.data;
        else
            error(me,'Projection item description missing');
        end
    end
    tag = find_tag(item,FIFF.FIFF_PROJ_ITEM_CH_NAME_LIST);
    if ~isempty(tag)
        namelist = tag.data;
    else
        error(me,'Projection item channel list missing');
    end
    tag = find_tag(item,FIFF.FIFF_PROJ_ITEM_KIND);
    if ~isempty(tag)
        kind = tag.data;
    else
        error(me,'Projection item kind missing');
    end
    tag = find_tag(item,FIFF.FIFF_PROJ_ITEM_NVEC);
    if ~isempty(tag)
        nvec = tag.data;
    else
        error(me,'Number of projection vectors not specified');
    end
    tag = find_tag(item,FIFF.FIFF_PROJ_ITEM_CH_NAME_LIST);
    if ~isempty(tag)
        names = fiff_split_name_list(tag.data);
    else
        error(me,'Projection item channel list missing');
    end
    tag = find_tag(item,FIFF.FIFF_PROJ_ITEM_VECTORS);
    if ~isempty(tag)
        data = tag.data;
    else
        error(me,'Projection item data missing');
    end
    tag = find_tag(item,FIFF.FIFF_MNE_PROJ_ITEM_ACTIVE);
    if ~isempty(tag)
        active = tag.data;
    else
        active = false;
    end
    if size(data,2) ~= length(names)
        error(me,'Number of channel names does not match the size of data matrix');
    end
    one.kind           = kind;
    one.active         = active;
    one.desc           = desc;
    %
    %   Use exactly the same fields in data as in a named matrix
    %
    one.data.nrow      = nvec;
    one.data.ncol      = nchan;
    one.data.row_names = [];
    one.data.col_names = names;
    one.data.data      = data;
    %
    projdata(i) = one;
end

if length(projdata) > 0
    fprintf(1,'\tRead a total of %d projection items:\n', ...
        length(projdata));
    for k = 1:length(projdata)
        fprintf(1,'\t\t%s (%d x %d)',projdata(k).desc, ...
            projdata(k).data.nrow,projdata(k).data.ncol);
        if projdata(k).active
            fprintf(1,' active\n');
        else
            fprintf(1,' idle\n');
        end
    end
end


return;

    function [tag] = find_tag(node,findkind)
        
        for p = 1:node.nent
            if node.dir(p).kind == findkind
                tag = fiff_read_tag(fid,node.dir(p).pos);
                return;
            end
        end
        tag = [];
    end


end
