function [lap,edge] = mesh_laplacian(vertex,face)
% MESH_LAPLACIAN: Laplacian of irregular triangular mesh
%
% Useage: [lap,edge] = mesh_laplacian(vertex,face)
%
% Returns 'lap', the Laplacian (2nd spatial derivative) of an
% irregular triangular mesh, and 'edge', the linear distances
% between vertices of 'face'.  'lap' and 'edge' are square,
% [Nvertices,Nvertices] in size, sparse in nature.
%
% It is assumed that 'vertex' contains the (x,y,z) Cartesian
% coordinates of each vertex and that 'face' contains the
% triangulation of vertex with indices into 'vertex' that
% are numbered from 1:Nvertices.  For information about
% triangulation, see 'help convhull' or 'help convhulln'.
%
% The neighbouring vertices of vertex 'i' is given by:
%
% k = find(edge(i,:));
%
% The math of this routine is given by:
%
% Oostendorp, Oosterom & Huiskamp (1989),
% Interpolation on a triangulated 3D surface.
% Journal of Computational Physics, 80: 331-343.
%
% See also, eeg_interp_scalp_mesh
%

% Licence:  GNU GPL, no implied or express warranties
% History:  04/2002, Darren.Weber@flinders.edu.au
%           - initial version was inefficient and incorrect
%             at one point.
%           (c) 04/2002 Robert Oostenveld
%           - completely revised/reconstructed code (lapcal.m)
%           - agreed to release into eeg_toolbox under GNU GPL
%           04/2002, Darren.Weber@flinders.edu.au
%           - modified edge initialization from sparse to
%             full matrix and slightly improved speed of
%             calculation for edge norms
%           - added tic/toc timing report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nvertex = size(vertex,1);
nface = size(face,1);

fprintf('MESH_LAPLACIAN: Calc Laplacian matrix for %5d vertices...',nvertex);
tic

% the matrix 'edge' is the connectivity of all vertices
edge = zeros(nvertex);
for i=1:nface,
   
    % compute the length of all triangle edges (Diff is [3x3])
    Diff = [vertex(face(i,[1 2 3]),:) - vertex(face(i,[2 3 1]),:)];
    Norm = sqrt( sum(Diff.^2, 2) );
   
    edge(face(i,1),face(i,2)) = Norm(1);
    edge(face(i,2),face(i,3)) = Norm(2);
    edge(face(i,3),face(i,1)) = Norm(3);
   
    % make sure that all edges are symmetric
    edge(face(i,2),face(i,1)) = Norm(1);
    edge(face(i,3),face(i,2)) = Norm(2);
    edge(face(i,1),face(i,3)) = Norm(3);
end

% Using edge to identify nearest vertices, calculate
% the Laplacian for an irregular mesh
lap = zeros(nvertex);
for i=1:nvertex,
   
    k = find(edge(i,:));        % the indices of the neighbours
   
    ni = length(k);             % the number of neighbours
   
    hi = mean(edge(i,k));       % the average distance to the neighbours
    invhi = mean(1./edge(i,k)); % the average inverse distance to the neighbours
   
    lap(i,i) = -(4/hi) * invhi; % Laplacian of vertex itself
   
    lap(i,k) =  (4/(hi*ni)) * 1./edge(i,k); % Laplacian of direct neighbours
   
    % Laplacian is zero for all indirect neighbours
    % See Oostendorp, Oosterom & Huiskamp (1989, pp. 334-335)
end

edge = sparse(edge);
lap = sparse(lap);

t = toc;
fprintf('done (%6.2f sec).\n',t);

return
