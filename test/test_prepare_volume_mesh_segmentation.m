function test_prepare_volume_mesh_segmentation

% test the function that generates meshes also used for constructing SIMBIO FEM head models
% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1815

% TEST test_prepare_volume_mesh
% TEST ft_datatype_segmentation 

%return;

%% creare an example segmentation 
example.dim = [91 104 111];   % slightly different numbers
example.transform = eye(4);   
example.coordsys = 'ctf';
example.unit = 'mm';
example.seg = zeros(example.dim);

% adjusting transformation matrix: center of the head-coordinates [0 0 0] should
% be the center of volume

% center of volume in voxel coordinates     
x = round(example.dim(1)/2);
y = round(example.dim(2)/2);
z = round(example.dim(3)/2);
    
x = round(x);
y = round(y);
z = round(z);      
             
origin = [x y z];

example.transform(1:4,4)   = [-origin(:); 1];  % head-coordinate [0 0 0] is in the center of 
                                                % the volume (x y z in voxel-coordinates)

% compute position for each voxel in voxelspace and in headspace

[X, Y, Z] = ndgrid(1:example.dim(1), 1:example.dim(2), 1:example.dim(3));
voxelpos = [X(:) Y(:) Z(:)];
headpos = warp_apply(example.transform, voxelpos);

% create 3 spheres

radius1 = 40;
radius2 = 30;
radius3 = 20;

for i=1:size(headpos,1)
    % from small to large
    if norm(headpos(i,:))<radius3
        example.seg(i) = 3;
    elseif norm(headpos(i,:))<radius2
        example.seg(i) = 2;
    elseif norm(headpos(i,:))<radius1
        example.seg(i) = 1;
    end
end

% plot it

% cfg=[];
% %cfg.interactive = 'yes';
% cfg.funparameter = 'seg';
% figure;
% ft_sourceplot(cfg,example);
%close all;

% convert it to probabilistic
seg=ft_datatype_segmentation(example,'segmentationstyle','probabilistic');


disp(seg)
%           dim: [91 104 111]
%     transform: [4x4 double]
%      coordsys: 'ctf'
%          unit: 'mm'
%      tissue_1: [91x104x111 logical]
%      tissue_2: [91x104x111 logical]
%      tissue_3: [91x104x111 logical]



%% test 
cfg=[];
%cfg.resolution = 1;
%cfg.tissue = 
mesh=prepare_volume_mesh_segmentation(cfg,seg);

 disp(mesh)
%             pnt: [283008x3 double]
%             hex: [267731x8 double]
%          tissue: [267731x1 double]
%     tissuelabel: {'101'  '102'  '103'}
%      tissuename: {'tissue_1'  'tissue_2'  'tissue_3'}

cfg=[];
%cfg.resolution = 1;
cfg.tissue = {'tissue_1' 'tissue_2'}; 
mesh2=prepare_volume_mesh_segmentation(cfg,seg);

disp(mesh2)
%             pnt: [253256x3 double]
%             hex: [234360x8 double]
%          tissue: [234360x1 double]
%     tissuelabel: {'101'  '102'}
%      tissuename: {'tissue_1'  'tissue_2'}

cfg=[];
%cfg.resolution = 1;
cfg.tissue = {'tissue_3'}; 
mesh3=prepare_volume_mesh_segmentation(cfg,seg);

  disp(mesh3)
%             pnt: [37224x3 double]
%             hex: [33371x8 double]
%          tissue: [33371x1 double]
%     tissuelabel: {'101'}
%      tissuename: {'tissue_3'}

cfg=[];
cfg.resolution = 2;
cfg.tissue = {'tissue_3'}; 
mesh4=prepare_volume_mesh_segmentation(cfg,seg);

disp(mesh4)
%             pnt: [37224x3 double]
%             hex: [166855x4 double]
%          tissue: [166855x1 double]
%     tissuelabel: {'1'}
%      tissuename: {'tissue_3'}

% vgrid is writing the wireframe mesh file, this may take some time ...
% WARNING: max_maxdim must divide nx,ny,nz!
% WARNING: enlarging image from size 91, 104, 111 to 92,106,112
% WARNING: gridGenerator: Changing grid generation to 'tetra5' 
% elapsed time: 3.231
% Trying to read FE mesh from file tp4f96190f_6895_4d28_b487_f938a19f7b37.v.
% Successfully read mesh with 37224 nodes and 166855 elements.
