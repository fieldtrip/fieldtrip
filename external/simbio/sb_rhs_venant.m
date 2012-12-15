function rhs = sb_rhs_venant(pos,dir,vol);

% SB_RHS_VENANT
%
% $Id$

%find node closest to source position
next_nd = sb_get_next_nd(pos,vol.pos);
%find nodes neighbouring closest node
if isfield(vol,'tet')
    ven_nd = sb_get_ven_nd(next_nd,vol.tet);
elseif isfield(vol,'hex')
    ven_nd = sb_get_ven_nd(next_nd,vol.hex);
else
    error('No connectivity information given!');
end
%calculate rhs matrix
loads = sb_calc_ven_loads(pos,dir,ven_nd,vol.pos);
%assign values in sparse matrix
i = reshape(ven_nd',[],1);
j = reshape(repmat([1:size(ven_nd,1)],size(ven_nd,2),1),[],1);
loads = reshape(loads',[],1);
j = j(i~=0);
loads = loads(i~=0);
i = i(i~=0);
rhs = sparse(i,j,loads,size(vol.pos,1),size(pos,1));
end

function next_nd = sb_get_next_nd(pos,node);
next_nd = zeros(size(pos,1),1);
for i=1:size(pos,1)
    [dist, next_nd(i)] = min(sum(bsxfun(@minus,node,pos(i,:)).^2,2));
end
end

function ven_nd = sb_get_ven_nd(next_nd,elem);
ven_nd = zeros(size(next_nd,1),1);
for i=1:size(next_nd,1)
    [tmp1,tmp2] = find(elem == next_nd(i));
    tmp = unique(elem(tmp1,:));
    %tmp = tmp(tmp~=next_nd(i)); seems like this is not done in the
    %original.
    if(length(tmp) > size(ven_nd,2))
        ven_nd = [ven_nd, zeros(size(ven_nd,1),length(tmp)-size(ven_nd,2))];
    end
    ven_nd(i,1:length(tmp)) = tmp;
end
end
