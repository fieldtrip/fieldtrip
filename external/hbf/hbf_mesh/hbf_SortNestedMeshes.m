function meshes_sorted=hbf_SortNestedMeshes(meshes)
%HBF_SORTNESTEDMESHES sorts nested triangle meshes ("small to large") 
%
% meshes_out=HBF_SORTNESTEDMESHES(meshes_in)
%   meshes_in:  the meshes to be sorted, cell array of hbf structs
%   meshes_out: sorted meshes, cell array of hbf structs
% 
% If the meshes are not nested, the function does not do anything.
% v160229 (c) Matti Stenroos

Nmesh=length(meshes);
prange=zeros(Nmesh,3);
for I=1:Nmesh
    prange(I,:)=max(meshes{I}.p)-min(meshes{I}.p);
end
orders=zeros(Nmesh,3);
for J=1:3
    [~,orders(:,J)]=sort(prange(:,J));
end
if ~all(orders(:,1)==orders(:,2)) && ~all(orders(:,1)==orders(:,3))
    fprintf('Meshes not nested. Not sorting...\n');
    meshes_sorted=meshes;
    return;
end
meshes_sorted=cell(Nmesh,1);
order=orders(:,1);
for I=1:Nmesh
    meshes_sorted{I}=meshes{order(I)};
end
if any (order-(1:Nmesh)')
    fprintf('Sorted meshes.\n');
else
    fprintf('Meshes were already in correct order.\n');   
end
