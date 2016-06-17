function [posR, triR, removetri] = remove_vertices(pos, tri, removepos)

% REMOVE_VERTICES removes specified vertices from a triangular mesh
% renumbering the vertex-indices for the triangles and removing all
% triangles with one of the specified vertices.
%
% Use as
%   [pos, tri] = remove_vertices(pos, tri, sel)

npos = size(pos,1);
ntri = size(tri,1);

if numel(removepos)==size(pos,1) && all(removepos==0 | removepos==1)
  removepos = find(removepos);
end

% remove the vertices and determine the new numbering (indices) in numb
keeppos = setdiff(1:npos, removepos);
numb    = zeros(1,npos);
numb(keeppos) = 1:length(keeppos);

% look for triangles referring to removed vertices
removetri = false(ntri,1);
removetri(ismember(tri(:,1), removepos)) = true;
removetri(ismember(tri(:,2), removepos)) = true;
removetri(ismember(tri(:,3), removepos)) = true;

% remove the vertices and triangles
posR = pos(keeppos, :);
triR = tri(~removetri,:);

% renumber the vertex indices for the triangles
triR = numb(triR);
