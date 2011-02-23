function [nodes] = fiff_dir_tree_find(tree, kind)
%
% [nodes] = fiff_dir_tree_find(tree,kind)
%
% Find nodes of the given kind from a directory tree structure
%

%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.4  2006/11/30 05:43:29  msh
%   Fixed help text in fiff_dir_tree_find
%   Fixed check for the existence of parent MRI block in mne_read_forward_solution
%
%   Revision 1.3  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.2  2006/04/12 10:29:02  msh
%   Made evoked data writing compatible with the structures returned in reading.
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%

me = 'MNE:fiff_dir_tree_find';

if nargin ~= 2
    error(me, 'Incorrect number of arguments');
end

nodes = struct('block', {}, 'id', {}, 'parent_id', {}, 'nent', {}, 'nchild', {}, 'dir', {}, 'children', {});
%
%   Am I desirable myself?
%
if tree.block == kind
    nodes = [ nodes tree ];
end
%
%   Search the subtrees
%
for k = 1:tree.nchild
    nodes = [ nodes fiff_dir_tree_find(tree.children(k), kind) ];
end

return

