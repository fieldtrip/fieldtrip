function fiff_list_dir_tree(out, tree, indent)

%
% fiff_list_dir_tree(fid, tree)
%
% List the fiff directory tree structure
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.3  2009/04/01 21:25:50  msh
%   Fixed problems with fiff nodes with no tags in them
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me = 'MNE:fiff_list_dir_tree';

if nargin == 2
    indent = 0;
end

for k = 1:indent
    fprintf(out, '\t');
end
if tree.block ~= 0
    fprintf(out, '{ %d\n', tree.block);
end

count = 1;
for k = 1:tree.nent
    if k == 1
        print = true;
        for p = 1:indent
            fprintf(out, '\t');
        end
        fprintf(out, 'tag : %d', tree.dir(k).kind);
    else
        if tree.dir(k).kind == tree.dir(k - 1).kind
            count = count + 1;
        else
            if count > 1
                fprintf(out, ' [%d]\n', count);
            else
                fprintf(out, '\n');
            end
            for p = 1:indent
                fprintf(out, '\t');
            end
            fprintf(out, 'tag : %d', tree.dir(k).kind);
            count = 1;
        end
    end
end
if count > 1
    fprintf(out, ' [%d]\n', count);
else
    fprintf(out, '\n');
end

for k = 1:tree.nchild
    fiff_list_dir_tree(out, tree.children(k), indent + 1);
end

for k = 1:indent
    fprintf(out, '\t');
end

if tree.block ~= 0
    fprintf(out, '} %d\n', tree.block);
end

return

