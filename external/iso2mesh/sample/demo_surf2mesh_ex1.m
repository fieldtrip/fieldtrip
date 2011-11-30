%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   demo script for mesh generation from surface patches and bounding box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preparation
% user must add the path of iso2mesh to matlab path list
% addpath('../');

% user need to add the full path to .../iso2mesh/bin directory
% to windows/Linux/Unix PATH environment variable

%% load the sample data
load tube_surface.mat

% f and v stores the surface patch faces and nodes
%% perform mesh generation
[node,elem,face]=surf2mesh(v,f,[1 1 1],[100 100 100],0.1,25);

%% visualize the resulting mesh
plotmesh(node,face(:,1:3));
axis equal;
