%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Create meshes for a sphere inside a cubic domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparation

% you have to add the path to iso2mesh toolbox 
% addpath('/path/to/iso2mesh/toolbox/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Part 0.  Create a Spherical Mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[no,el]=meshasphere([30 30 30],20,2.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Part I.  A Coarse Mesh for a Sphere Inside a Box with Refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate a coarse volumetric mesh from the sphere with an additional bounding box
% the maximum element volume is 8

ISO2MESH_SESSION='demo_sph3_';

srcpos=[30. 30. 0.];                     % set the center of the ROI
fixednodes=[30.,30.,0.1; 30 30 30];     % add control points so we can refine mesh densities there
nodesize=[ones(size(no,1),1) ; 0.2; 4];  % set target edge size of 1 for all nodes on the sphere
                                         % target edge size 0.3 near (30,30,0.05)
                                         % and target edge size 4 near (30,30,30)
nfull=[no;fixednodes];                   % append additional control points
[node3,elem3,face3]=surf2mesh([nfull,nodesize],el,[0 0 0],[61 61 61],1,8,[30 30 30],[],[2 2 2 2 6 6 6 6]);
                             % ^- add node size as the last            ^ max volume     ^- edge sizes at the 8 
                             %    column to node                                           corners of the bounding box
[node3,elem3,face3]=sortmesh(srcpos,node3,elem3,1:4,face3,1:3);  % reorder the nodes/elements 
                                                 % so that the nodes near earch order
                                                 % are more clustered in the memory
elem3(:,1:4)=meshreorient(node3,elem3(:,1:4));   % reorient elements to ensure the volumns are positive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Part II.  A Dense Mesh for a Sphere Inside a Box with Refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a dense volumetric mesh from the sphere with an additional bounding box
% the maximum element volume is 2

ISO2MESH_SESSION='demo_sph2_';

nodesize=[0.7*ones(size(no,1),1) ; 0.2; 2];  % set target edge size to 0.7 near the sphere
                                             % 0.2 near (30,30,0.5) and 2 near (30,30,30)
[node2,elem2,face2]=surf2mesh([nfull,nodesize],el,[0 0 0],[61 61 61],1,2,[30 30 30],[],[1 1 1 1 5 5 5 5]);

figure; plotmesh(node2,face2(:,1:3),'y>30');axis equal;

[node2,elem2,face2]=sortmesh(srcpos,node2,elem2,1:4,face2,1:3);
elem2(:,1:4)=meshreorient(node2,elem2(:,1:4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Part III.  A Coarse Mesh for a Sphere Inside a Box without Refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ISO2MESH_SESSION='demo_sph1_';

% reduce the surface node numbers to 20%  
[no2,el2]=meshresample(no,el,0.2);  % down sample the sphere mesh

% using the coarse spherical surface, we generate a coarse volumetric
% mesh with maximum volume of 10

[node1,elem1,face1]=surf2mesh(no2,el2,[0 0 0],[61 61 61],1,10,[30 30 30],[],1);
[node1,elem1,face1]=sortmesh(srcpos,node1,elem1,1:4,face1,1:3);
elem1(:,1:4)=meshreorient(node1,elem1(:,1:4));

clear ISO2MESH_SESSION
