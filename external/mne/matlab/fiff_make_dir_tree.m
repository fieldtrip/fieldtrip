function [tree, last] = fiff_make_dir_tree(fid,dir,start,indent)
%
% [tree, last] = fiff_make_dir_tree(fid,dir,start,indent)
%
% Create the directory tree structure
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.6  2009/04/01 21:25:50  msh
%   Fixed problems with fiff nodes with no tags in them
%
%   Revision 1.5  2006/04/26 19:50:58  msh
%   Added fiff_read_mri
%
%   Revision 1.4  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
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
%  Define the relevant constants here
%  no need to get the whole fiff_define_costants
%
FIFF_BLOCK_START     = 104;
FIFF_BLOCK_END       = 105;
FIFF_FILE_ID         = 100;
FIFF_BLOCK_ID        = 103;
FIFF_PARENT_BLOCK_ID = 110;

me='MNE:fiff_make_dir_tree';

verbose=0;

if nargin == 2
    indent = 0;
    start = 1;
elseif nargin == 3
    indent = 0;
elseif nargin ~= 4
    error(me,'Incorrect number of arguments');
end

if dir(start).kind == FIFF_BLOCK_START
    tag = fiff_read_tag(fid,dir(start).pos);
    block = tag.data;
else
    block = 0;
end

if verbose ~= 0
    for k = 1:indent
        fprintf(1,'\t');
    end
    fprintf(1,'start { %d\n',block);
end

nchild = 0;
this = start;

tree.block    = block;
tree.id        = [];
tree.parent_id = [];
tree.nent     = 0;
tree.nchild   = 0;
tree.dir      = dir(this);
tree.children = struct('block', {}, 'id', {}, 'parent_id', {}, 'nent', {}, 'nchild', {}, 'dir', {}, 'children', {});
while this <= length(dir)
    if dir(this).kind == FIFF_BLOCK_START
        if this ~= start
            [ child , this ] = fiff_make_dir_tree(fid,dir,this,indent+1);
            tree.nchild = tree.nchild + 1;
            tree.children(tree.nchild) = child;

        end
    elseif dir(this).kind == FIFF_BLOCK_END
        tag = fiff_read_tag(fid,dir(start).pos);
        if tag.data == block
            break;
        end
    else
        tree.nent = tree.nent + 1;
        tree.dir(tree.nent) = dir(this);
        %
        %  Add the id information if available
        %
        if block == 0
            if dir(this).kind == FIFF_FILE_ID
                tag = fiff_read_tag(fid,dir(this).pos);
                tree.id = tag.data;
            end
        else
            if dir(this).kind == FIFF_BLOCK_ID
                tag = fiff_read_tag(fid,dir(this).pos);
                tree.id = tag.data;
            elseif dir(this).kind == FIFF_PARENT_BLOCK_ID
                tag = fiff_read_tag(fid,dir(this).pos);
                tree.parent_id = tag.data;
            end
        end
    end
    this = this + 1;
end
%
% Eliminate the empty directory
%
if tree.nent == 0
    tree.dir = [];
end
if verbose ~= 0
    for k = 1:indent+1
        fprintf(1,'\t');
    end
    % fprintf(1,'block = %d nent = %d nchild = %d\n',tree.block,tree.nent,tree.nchild);
    fprintf(1,'block=%d\n',tree.block);
    for k = 1:indent
        fprintf(1,'\t');
    end
    fprintf(1,'end } %d\n',block);
end

last = this;

return;

