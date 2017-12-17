function test_ft_plot_mesh

% MEM 1gb
% WALLTIME 00:15:00

% TEST test_ft_plot_mesh
% TEST ft_plot_mesh

%% the first is a simple triangluar mesh

bnd.pos = randn(20,3);
bnd.tri = delaunay(bnd.pos(:,1:2));

figure
ft_plot_mesh(bnd);

%% the following was made by Lilla Magyari in 2013 and pertains to FEM meshes

% create an example segmentation
example.dim = [91 104 111]; % slightly different numbers
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

example.transform(1:4,4) = [-origin(:); 1]; % head-coordinate [0 0 0] is in the center of
% the volume (x y z in voxel-coordinates)

% compute position for each voxel in voxelspace and in headspace

[X, Y, Z] = ndgrid(1:example.dim(1), 1:example.dim(2), 1:example.dim(3));
voxelpos = [X(:) Y(:) Z(:)];
headpos = ft_warp_apply(example.transform, voxelpos);

% create 3 spheres

radius1 = 40;
radius2 = 30;
radius3 = 20;

for i = 1:size(headpos,1)
  % from small to large
  if norm(headpos(i,:))<radius3
    example.seg(i) = 3;
  elseif norm(headpos(i,:))<radius2
    example.seg(i) = 2;
  elseif norm(headpos(i,:))<radius1
    example.seg(i) = 1;
  end
end


% convert it to probabilistic
seg = ft_datatype_segmentation(example,'segmentationstyle','probabilistic');

% create smaller segmentation

example.dim = [10 11 9]; % slightly different numbers
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

example.transform(1:4,4) = [-origin(:); 1]; % head-coordinate [0 0 0] is in the center of
% the volume (x y z in voxel-coordinates)

% compute position for each voxel in voxelspace and in headspace

[X, Y, Z] = ndgrid(1:example.dim(1), 1:example.dim(2), 1:example.dim(3));
voxelpos = [X(:) Y(:) Z(:)];
headpos = ft_warp_apply(example.transform, voxelpos);

% create 3 spheres

radius1 = 4;
radius2 = 3;
radius3 = 2;

for i = 1:size(headpos,1)
  % from small to large
  if norm(headpos(i,:))<radius3
    example.seg(i) = 3;
  elseif norm(headpos(i,:))<radius2
    example.seg(i) = 2;
  elseif norm(headpos(i,:))<radius1
    example.seg(i) = 1;
  end
end

% convert it to probabilistic
seg_small = ft_datatype_segmentation(example,'segmentationstyle','probabilistic');

%% create mesh
cfg = [];
cfg.method = 'hexahedral';
%cfg.tissue =
mesh = ft_prepare_mesh(cfg,seg);

% 1 cube mesh
mesh0.pos = mesh.pos;
mesh0.hex = mesh.hex(1,:);
mesh0.tissue = mesh.tissue(1,:);
mesh0.tissuelabel = mesh.tissuelabel;
mesh0.unit = mesh.unit;

% few hundreds cube mesh

cfg = [];
cfg.method = 'hexahedral';
%cfg.tissue =
mesh2 = ft_prepare_mesh(cfg,seg_small);

%% plot mesh

% larger mesh
figure; ft_plot_mesh(mesh,'surfaceonly','yes');
% smaller mesh
figure; ft_plot_mesh(mesh2,'surfaceonly','yes');
% smallest mesh
figure; ft_plot_mesh(mesh0,'surfaceonly','yes');

