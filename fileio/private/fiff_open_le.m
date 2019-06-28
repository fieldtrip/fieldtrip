function [fid, tree, dir] = fiff_open_le(fname)
%
% [fid, tree, dir] = fiff_open_le(fname)
%
% Open a fif file and provide the directory of tags
%
% fid     the opened file id
% tree    tag directory organized into a tree
% dir     the sequential tag directory
%
% This is a modified version, specific for opening 'little endian' fiff files
% Arjen Stolk

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.7  2009/03/30 11:37:37  msh
%   Added copying of measurement info blocks from the original like in mne_browse_raw
%
%   Revision 1.6  2008/11/16 21:31:23  msh
%   Added mne_transform_coordinates and new coordinate frame definitions
%
%   Revision 1.5  2006/05/03 19:03:19  msh
%   Eliminated the use of cast function for Matlab 6.5 compatibility
%
%   Revision 1.4  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.3  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.2  2006/04/17 15:01:34  msh
%   More small improvements.
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_open';
verbose=false;
fid = fopen_or_error(fname,'rb','ieee-le'); % Arjen Stolk: this is 'ieee-be' in fiff_open.m

%
%   Check that this looks like a fif file
%
tag = fiff_read_tag_info(fid);
if tag.kind ~= FIFF.FIFF_FILE_ID
    ft_error(me,'file does not start with a file id tag');
end
if tag.type ~= FIFF.FIFFT_ID_STRUCT
    ft_error(me,'file does not start with a file id tag');
end
if tag.size ~= 20
    ft_error(me,'file does not start with a file id tag');
end
tag = fiff_read_tag(fid);
if tag.kind ~= FIFF.FIFF_DIR_POINTER
    ft_error(me,'file does have a directory pointer');
end
if nargout == 1
    fseek(fid,0,'bof');
    return;
end
%
%   Read or create the directory tree
%
if verbose
    fprintf(1,'\tCreating tag directory for %s...',fname);
end
dirpos = double(tag.data);
if dirpos > 0
    tag = fiff_read_tag(fid,dirpos);
    dir = tag.data;
else
    k = 0;
    fseek(fid,0,'bof');
    dir = struct('kind',{},'type',{},'size',{},'pos',{});
    while tag.next >= 0
        pos = ftell(fid);
        tag = fiff_read_tag_info(fid);
        k = k + 1;
        dir(k).kind = tag.kind;
        dir(k).type = tag.type;
        dir(k).size = tag.size;
        dir(k).pos  = pos;
    end
end
%
%   Create the directory tree structure
%
tree = fiff_make_dir_tree(fid,dir);
if verbose
    fprintf(1,'[done]\n');
end
%
%   Back to the beginning
%
fseek(fid,0,'bof');
return;
