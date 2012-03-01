function [edges,idx,edgemap]=uniqedges(elem)

if(size(elem)==2)
   edges=elem;
elseif(size(elem)>=3)
   edges=meshedge(elem);
else
   error('invalid input');
end

[uedges,idx,jdx]=unique(sort(edges,2),'rows');
edges=edges(idx,:);
edgemap=reshape(jdx,size(elem));
