function fiff_copy_tree(fidin, in_id, nodes, fidout)
%
%    fiff_copy_tree(fidin, in_id, nodes, fidout)
%
%    Copies directory subtrees from fidin to fidout
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%    Revision 1.2  2009/04/01 21:25:50  msh
%    Fixed problems with fiff nodes with no tags in them
%
%    Revision 1.1  2009/03/30 11:37:37  msh
%    Added copying of measurement info blocks from the original like in mne_browse_raw
%
%

me = 'MNE:fiff_copy_tree';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin ~= 4
    error(me, 'Incorrect number of arguments');
end

if length(nodes) <= 0
    return
end

FIFFV_NEXT_SEQ = 0;

for k = 1:length(nodes)
    fiff_start_block(fidout, nodes(k).block);
    if ~isempty(nodes(k).id)
        if ~isempty(in_id)
            fiff_write_id(fidout, FIFF.FIFF_PARENT_FILE_ID, in_id);
        end
        fiff_write_id(fidout, FIFF.FIFF_BLOCK_ID);
        fiff_write_id(fidout, FIFF.FIFF_PARENT_BLOCK_ID, nodes(k).id);
    end
    for p = 1:nodes(k).nent
        %
        %   Do not copy these tags
        %
        if nodes(k).dir(p).kind == FIFF.FIFF_BLOCK_ID || nodes(k).dir(p).kind == FIFF.FIFF_PARENT_BLOCK_ID || nodes(k).dir(p).kind == FIFF.FIFF_PARENT_FILE_ID
            continue;
        end
        %
        %   Read and write tags, pass data through transparently
        %
        if fseek(fidin, nodes(k).dir(p).pos, 'bof') == -1
            error(me, 'Could not seek to the tag');
        end
        tag.kind = fread(fidin, 1, 'int32');
        tag.type = fread(fidin, 1, 'uint32');
        tag.size = fread(fidin, 1, 'int32');
        tag.next = fread(fidin, 1, 'int32');
        tag.data = fread(fidin, tag.size, 'uchar');
        count = fwrite(fidout, tag.kind, 'int32');
        if count ~= 1
            error(me, 'write failed.');
        end
        count = fwrite(fidout, tag.type, 'int32');
        if count ~= 1
            error(me, 'write failed.');
        end
        count = fwrite(fidout, tag.size, 'int32');
        if count ~= 1
            error(me, 'write failed');
        end
        count = fwrite(fidout, int32(FIFFV_NEXT_SEQ), 'int32');
        if count ~= 1
            error(me, 'write failed');
        end
        count = fwrite(fidout, tag.data, 'uchar');
        if count ~= tag.size
            error(me, 'write failed');
        end
    end
    for p = 1:nodes(k).nchild
        fiff_copy_tree(fidin, in_id, nodes(k).children(p), fidout);
    end
    fiff_end_block(fidout, nodes(k).block);
end

return


