function [edges] = mne_mesh_edges(faces)

% MESH_EDGES   Returns sparse matrix with edges number
%
%   SYNTAX
%       [EDGES] = MESH_EDGES(FACES)
%
%   faces : matrix of size [n_trianges, 3]
%   edges : sparse matrix of size [n_vertices, n_vertices]

%
%   Author : Alexandre Gramfort, MGH Martinos Center
%   License : BSD 3-clause
%

me = 'MESH_EDGES';

if nargin == 0
    eval(['help ',lower(me)])
    return
end

faces = double(faces);
npoints = max(faces(:));
nfaces = size(faces,1);
edges = sparse(faces(:,1),faces(:,2),ones(1,nfaces),npoints,npoints);
edges = edges + sparse(faces(:,2),faces(:,3),ones(1,nfaces),npoints,npoints);
edges = edges + sparse(faces(:,3),faces(:,1),ones(1,nfaces),npoints,npoints);

edges = edges + edges';

end %  function
