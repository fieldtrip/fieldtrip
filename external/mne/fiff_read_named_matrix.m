function [ mat ] = fiff_read_named_matrix(fid,node,matkind)

%
% [mat] = fiff_read_named_matrix(fid,node)
%
% Read named matrix from the given node
%


%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.4  2007/11/13 10:55:32  msh
%   Specify the second argument to all calls to the exist function
%
%   Revision 1.3  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.2  2006/04/20 21:49:38  msh
%   Added mne_read_inverse_operator
%   Changed some of the routines accordingly for more flexibility.
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_read_named_matrix';

if nargin ~= 3
    error(me,'Incorrect number of arguments');
end
%
%   Descend one level if necessary
%
found_it=false;
if node.block ~= FIFF.FIFFB_MNE_NAMED_MATRIX
    for k = 1:node.nchild
        if node.children(k).block == FIFF.FIFFB_MNE_NAMED_MATRIX
            if has_tag(node.children(k),matkind)
                node = node.children(k);
                found_it = true;
                break;
            end
        end
    end
    if ~found_it
        error(me,'Desired named matrix (kind = %d) not available',matkind);
    end
else
    if ~has_tag(node,matkind)
        error(me,'Desired named matrix (kind = %d) not available',matkind);
    end
end
%
%   Read everything we need
%
tag = find_tag(node,matkind);
if isempty(tag)
    error(me,'Matrix data missing');
else
    data = tag.data;
end
nrow = size(data,1);
ncol = size(data,2);
tag = find_tag(node,FIFF.FIFF_MNE_NROW);
if ~isempty(tag)
    if tag.data ~= nrow
        error(me,'Number of rows in matrix data and FIFF_MNE_NROW tag do not match');
    end
end
tag = find_tag(node,FIFF.FIFF_MNE_NCOL);
if ~isempty(tag)
    if tag.data ~= ncol
        error(me,'Number of columns in matrix data and FIFF_MNE_NCOL tag do not match');
    end
end
tag = find_tag(node,FIFF.FIFF_MNE_ROW_NAMES);
if ~isempty(tag)
    row_names = tag.data;
end
tag = find_tag(node,FIFF.FIFF_MNE_COL_NAMES);
if ~isempty(tag)
    col_names = tag.data;
end
%
%   Put it together
%
mat.nrow = nrow;
mat.ncol = ncol;
if exist('row_names','var')
    mat.row_names = fiff_split_name_list(row_names);
else
    mat.row_names = [];
end
if exist('col_names','var')
    mat.col_names = fiff_split_name_list(col_names);
else
    mat.col_names = [];
end
mat.data = data;

return;


    function [tag] = find_tag(node,findkind)
        
        for p = 1:node.nent
            if node.dir(p).kind == findkind
                tag = fiff_read_tag(fid,node.dir(p).pos);
                return;
            end
        end
        tag = [];
        return;
    end

    function [has] = has_tag(this,findkind)
        
        for p = 1:this.nent
            if this.dir(p).kind == findkind
                has = true;
                return;
            end
        end
        has = false;
        return;
        
    end

end
