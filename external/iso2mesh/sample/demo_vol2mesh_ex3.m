% sample script to create volumetric mesh from 
% multiple levelsets of a binary segmented head image.
%
% Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>

% load full head image (T1 MRI scan)

fprintf(1,'loading binary head image...\n');
for i=1:256 
  head(:,:,i)=imread('head.tif',i);
end

% load segmented brain images (by freesurfer recon_all)
fprintf(1,'loading binary brain image...\n');
for i=1:256
  brain(:,:,i)=imread('brain.tif',i);
end

% fill holes in the head image and create the canonical binary volume
% this may take a few minutes for a 256x256x256 volume
fprintf(1,'filling holes in the volumetric images...\n');
tic
cleanimg=fillholes3d(logical(head>0),10);
cleanbrain=fillholes3d(logical(brain>0),10);
toc

% add brain image as additional segment
cleanimgfull=cleanimg+(cleanbrain>0);

% create volumetric tetrahedral mesh from the two-layer 3D images
% this may take another few minutes for a 256x256x256 volume
clear opt;

% set method for vol2mesh to 'simplify' to use these option
opt(1).keepratio=0.05; % resample levelset 1 to 5%  
opt(2).keepratio=0.1;  % resample levelset 2 to 10%

% by default, vol2mesh uses 'cgalsurf' method, which requires the following
opt(1).radbound=4; % head surface element size bound
opt(2).radbound=2; % brain surface element size bound
opt(1).side='lower'; %
opt(2).side='lower'; %

tic
[node,elem,face]=vol2mesh(cleanimgfull,1:size(cleanimg,1),1:size(cleanimg,2),1:size(cleanimg,3),opt,100,1);
toc

% plot the boundary surface of the generated mesh
h=slice(cleanimgfull,[],[120],[120 180]);
set(h,'linestyle','none')
hold on

plotmesh(node(:,[2 1 3]),face,'facealpha',0.7);
