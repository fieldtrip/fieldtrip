%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hello World - The "getting started" example of iso2mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this example, we illustrate the basic mesh density and precision 
% controls in iso2mesh.

%% the first step is to load a binary or gray-scale image 

hw=imread('helloworld.png');

% now we stack the image and form a binary volume

hw3d=1-repmat(hw,[1 1 50]);

%% create 3D mesh          |--------> threshold at v=0.7
[node,elem,face]=v2m(hw3d,0.7,5,40);
%                             |  |-> maximum volume
%                             |----> maximum surface triangle size

figure;
subplot(211);
plotmesh(node,face);axis equal;view(90,60);
subplot(212);
plotmesh(node,elem,'z<20');axis equal;view(90,60);

% mesh with denser surface    |----> surface triangle size is now 2
[node,elem,face]=v2m(hw3d,0.7,2,40);

figure;
subplot(211);
plotmesh(node,face);axis equal;view(90,60);
subplot(212);
plotmesh(node,elem,'z<20');axis equal;view(90,60);

%% create 3D mesh from gray-scale image to get smoother boundary
hw=imread('helloworld_gray.png');
hw3d=255-repmat(hw,[1 1 50]);
[node,elem,face]=v2m(hw3d,128,2,40);

figure;
subplot(211);
plotmesh(node,face);axis equal;view(90,60);
subplot(212);
plotmesh(node,elem,'z<20');axis equal;view(90,60);

%% create 3D mesh from gray-scale image with advanced distbound control
clear opt
opt.radbound=4;       % set surface triangle maximum size
opt.distbound=0.2;    % set max distance that deviates from the level-set
opt.autoregion=1;
[node,elem,face]=v2m(hw3d,128,opt,40);

figure;
subplot(211);
plotmesh(node,face);axis equal;view(90,60);
subplot(212);
plotmesh(node,elem,'z<20');axis equal;view(90,60);

