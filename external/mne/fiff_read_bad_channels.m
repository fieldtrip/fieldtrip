function [bads] = fiff_read_bad_channels(fid,node,ch_rename)
%
% [bads] = fiff_read_bad_channels(fid,node[,ch_rename])
%
% Reas the bad channel list from a node if it exists
%
% fid       - The file id
% node      - The node of interes
% ch_rename - Short-to-long channel name mapping
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.5  2006/05/03 18:53:04  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.4  2006/05/03 15:51:17  msh
%   Re-added fiff_read_bad_channels
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/12 17:09:28  msh
%   Added routines for reading noise-covariance matrices
%
%

me='MNE:fiff_read_bad_channels';

if nargin == 2
    ch_rename = {};
elseif nargin ~= 3
    error(me,'Incorrect number of arguments');
end

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

node = fiff_dir_tree_find(node,FIFF.FIFFB_MNE_BAD_CHANNELS);

bads = [];
if ~isempty(node)
    tag = find_tag(node,FIFF.FIFF_MNE_CH_NAME_LIST);
    if ~isempty(tag)
        bads = fiff_split_name_list(tag.data);
    end
end
bads = fiff_rename_list(bads, ch_rename);

return;

    function [tag] = find_tag(node,findkind)

        for p = 1:node.nent
            kind = node.dir(p).kind;
            pos  = node.dir(p).pos;
            if kind == findkind
                tag = fiff_read_tag(fid,pos);
                return;
            end
        end
        tag = [];
        return
    end

end
