function [stiff diinsy cols sysmat] = sb_calc_stiff(node,elem,cond);
if(~(size(node,2)==3))
    if(size(node,1)==3)
        node = node';
        warning('Dimension of wf.vol.nd should be #nodes x 3!')
    else
        error('wf.vol.nd has wrong dimension!')
    end
end
nn = size(node,1);
nn = int32(nn);
%node = [node(:,1);node(:,2);node(:,3)];
if(min(elem(:))==0)
    elem = elem + 1;
    warning('Numbering of nodes must start at 1 (Fortran numbering)!')
elseif (~min(min(elem))<0)
    error('No negative indices for cell definition allowed!')
end
if((size(elem,1)==4)|(size(elem,1)==8))
    mele = size(elem,1);
elseif((size(elem,2)==4)|(size(elem,2)==8))
    elem = elem';
    mele = size(elem,1);
    warning('Dimension of wf.vol.el should be (#nodes per ele) x #elem!')
else
    error('wf.vol.el has wrong dimension!')
end
if(mele==4)
    elem = [elem; zeros(4,size(elem,2))];
elseif(~(mele==8))
    error('Invalid number of nodes per element!')
end
if(~((size(cond,2)==1|size(cond,2)==6)&(size(cond,1)==size(elem,2))))
    if((size(cond,1)==1|size(cond,1)==6)&(size(cond,2)==size(elem,2)))
        cond = cond';
        warning('Dimension of wf.vol.cond should be #elem x #aniso!')
    else
        error('wf.vol.cond has wrong dimension!')
    end
end
%check für cond einfügen!
mele = int32(mele);
elem = int32(elem);
[diinsy,cols,sysmat] = calc_stiff_matrix_val(node,elem,cond,mele);
nn = double(nn);
diinsy = double(diinsy);
cols = double(cols);
rows = sb_sparse_to_mat(diinsy);
stiff = sparse(rows,cols,sysmat,nn,nn,length(sysmat));
end