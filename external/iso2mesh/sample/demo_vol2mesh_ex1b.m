%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   demo script for mesh generation from binarized volumetric image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preparation
% user must add the path of iso2mesh to matlab path list
% addpath('../');

% user need to add the full path to .../iso2mesh/bin directory
% to windows/Linux/Unix PATH environment variable

%% load the sample data
load rat_head.mat

% volimage is a volumetric image such as an X-ray or MRI image
% A,b are registration matrix and vector, respectively
%% perform mesh generation

%% use the alternative 'simplify' method: first create voxel-based
% surface mesh, and then resample it to desired density.
% this method does not guarantee to be free of self-intersecting
% element, as 'cgalsurf' promises.

[node,elem,face]=vol2mesh(volimage>0.05,1:size(volimage,1),1:size(volimage,2),...
                          1:size(volimage,3),0.1,2,1,'simplify');

%% visualize the resulting mesh

plotmesh(node,face);
axis equal;
