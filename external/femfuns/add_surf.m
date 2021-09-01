function newsurf = add_surf(surface1, surface2)


% ADD_SURF combine two surfaces and remove duplicate nodes and faces
% surface1 is the base, new nodes and faces of surface2 are added to it
% NOTE: self-intersections or non-connected surfaces are not resolved

%add new nodes
flag = ~ismember(surface2.node,surface1.node,'rows');
newsurf.node = [surface1.node; surface2.node(flag,:)];
newsurf.face = surface1.face;

%loop over faces of the to be added surface
for i = 1:size(surface2.face,1)
    %find the new face indices
    tmpfc = [];
    for j = 1:3
        [flag,idx] = ismember(surface2.node(surface2.face(i,j),:), newsurf.node, 'rows');
        if length(idx) == 1
            tmpfc = [tmpfc idx];
        else
            error('Only one vertex coordinate index should be found. Remove duplicates first')
        end
    end

    %check if it is a new face, then add
    flag = ~ismember(sort(tmpfc),sort(newsurf.face,2),'rows');
    if sum(flag)==1
        newsurf.face = [newsurf.face;tmpfc];
    end
end
