function [stiff diinsy cols sysmat] = sb_calc_stiff(vol)

% SB_CALC_STIFF
%
% $Id$

if(~(size(vol.pos,2)==3))
    if(size(vol.pos,1)==3)
        node = vol.pos';
        warning('Dimension of vol.pos should be #nodes x 3!')
    else
        error('vol.pos has wrong dimension!')
    end
else
    node = vol.pos;
end
npnt = size(node,1);
npnt = int32(npnt);

if isfield(vol,'tet')
    if size(vol.tet,1) == 4
        mele = size(vol.tet,1);
        elem = vol.tet;
    elseif size(vol.tet,2) == 4
        mele = size(vol.tet,2);
        elem = vol.tet';
    else
        error('vol.tet has wrong dimensions!')
    end
    elem = [elem; zeros(4,size(elem,2))];
elseif isfield(vol,'hex')
    if size(vol.hex,1) == 8
        mele = size(vol.hex,1);
        elem = vol.hex;
    elseif size(vol.hex,2) == 8
        mele = size(vol.hex,2);
        elem = vol.hex';
    else
        error('vol.hex has wrong dimensions!')
    end
else
    error('Could not find connectivity information!')
end

if min(min(elem(1:mele,:))) == 0
    elem = elem + 1;
    warning('Numbering of nodes in vol.tet/vol.hex must start at 1 (Fortran numbering)!')
elseif min(min(elem(1:mele,:))) < 0
    error('No negative indices for conectivity information allowed!')
end
    
if isfield(vol,'cond') && isfield(vol,'tissue') && isfield(vol,'tissuelabel')
    if length(vol.tissuelabel) == length(vol.cond)
         if length(vol.tissue) == size(elem,2)
             cond = zeros(size(elem,2),1);
             numlabels = length(vol.tissuelabel);
             for i=1:numlabels
                 cond(vol.tissue == i) = vol.cond(i);
             end
        else
            error('Dimensions of vol.tet or vol.hex and vol.tissue do not fit!');
        end
    else
        error('Dimensions of vol.cond and entries of vol.tissuelabel do not fit!');
    end
end

mele = int32(mele);
elem = int32(elem);

% check whether the nodes have right orientation

if isfield(vol,'tet')
    if ~sb_test_ori(node,elem(1:4,:)')
        error('Elements have wrong orientation, consider exchanging node 3 and 4');
        return;
    end
elseif isfield(vol,'hex')
    if ~sb_test_ori(node,elem')
        error('Elements have wrong orientation or are degenerated');
        return
    end
end

try
    [diinsy,cols,sysmat] = calc_stiff_matrix_val(node,elem,cond,mele);
catch err
    if ispc && strcmp(err.identifier,'MATLAB:invalidMEXFile')
        error('Error executing mex-file. Microsoft Visual C++ 2008 Redistributables and Intel Visual Fortran Redistributables are required.')
    else
        rethrow(err)
    end
end
npnt = double(npnt);
diinsy = double(diinsy);
cols = double(cols);
rows = sb_sparse_to_mat(diinsy);
stiff = sparse(rows,cols,sysmat,npnt,npnt,length(sysmat));
end
