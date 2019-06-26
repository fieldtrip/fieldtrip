function [pos] = moveinward(pos, move)
%This functions moves 'pos' inward according to their normals by 'move'
%units
propos = elproj(pos); % projection to 2D
tri = delaunay(propos); %creates delaunay triangulation of 2D plane, which will be used for the the 3D case
nor = normals(pos,tri); %compute normals of surface
ori = surfaceorientation(pos, tri, nor);
if ori==1
  % the normals are outward oriented
elseif ori==-1
  % the normals are inward oriented
  nor = -nor;
else
  ft_warning('cannot determine the orientation of the vertex normals');
end
pos = pos-move*nor; % moves pos inwards according to their normals
end

