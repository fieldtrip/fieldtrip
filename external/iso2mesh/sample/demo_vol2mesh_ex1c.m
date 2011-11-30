%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   demo script for mesh generation from binary volumetric image
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

%% use the alternative 'cgalmesh' method. This will call 
% cgalmesher to process labled volume to produce surfaces
% and tetrahedral mesh in a single run.
clear opt
opt.radbound=2;
[node,elem,face]=v2m(uint8(volimage),0.5,opt,100,'cgalmesh');


%% visualize the resulting mesh

plotmesh(node,face(:,1:3));
axis equal;
