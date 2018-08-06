function test_bug1836

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_datatype_segmentation ft_prepare_mesh ft_write_headshape

% this bugscript does not need to run automatically, because the problem is
% that ft_write_headshape does not write hexahedral meshes in the proper
% way into a ply format file. This can be detected only when the ply file
% is read into paraview.

%% make segmentation with a small size
% when the segmentation is small, the mesh will contain a few elements
% which is better for visualization

% example segmentation
example.dim = [91 104 111];   % slightly different numbers
example.transform = eye(4);
example.coordsys = 'ctf';
example.unit = 'mm';
example.seg = zeros(example.dim);

x = round(example.dim(1)/2);
y = round(example.dim(2)/2);
z = round(example.dim(3)/2);

x = round(x);
y = round(y);
z = round(z);

origin = [x y z];

example.transform(1:4,4)   = [-origin(:); 1];

[X, Y, Z] = ndgrid(1:example.dim(1), 1:example.dim(2), 1:example.dim(3));
voxelpos = [X(:) Y(:) Z(:)];
headpos = ft_warp_apply(example.transform, voxelpos);

% create a sphere with 1 and 2 mm radius in the example volume

example1 = example;
example2 = example;
radius1 = 1;
radius2 = 2;

for i=1:size(headpos,1)
  % from small to large
  if norm(headpos(i,:))<radius1
    example1.seg(i) = 1;
  end
  if norm(headpos(i,:))<radius2
    example2.seg(i) = 1;
  end
  
end
clear('example','origin','voxelpos','headpos','X','Y','Z','x','y','z','radius1','radius2');

% convert it to probabilistic
example1=ft_datatype_segmentation(example1,'segmentationstyle','probabilistic');
example2=ft_datatype_segmentation(example2,'segmentationstyle','probabilistic');


%% make meshes
cfg=[];
cfg.method = 'hexahedral';
%cfg.tissue =
mesh1=ft_prepare_mesh(cfg,example1);
mesh2=ft_prepare_mesh(cfg,example2);

% plot mesh in FieldTrip
figure
ft_plot_mesh(mesh1,'facelalpha',0.5); % single cube, rotate it to see it well
figure;
ft_plot_mesh(mesh2,'facealpha',0.5); % rubiks cube, rotate it to see it well

%% write meshes
% change path
ft_write_headshape('mesh1cube',mesh1,'format','ply')
ft_write_headshape('mesh27cube',mesh2,'format','ply')

%% do the same with a tetraheder
mesh = [];
mesh.pnt = [
  0 0 0
  1 0 0
  0 1 0
  0 0 1];
mesh.tet = [1 2 3 4];

ft_write_headshape('mesh1tet',mesh,'format','ply')

%% open it with paraview or meshlab and compare it to FT plot



