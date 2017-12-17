function test_bug3196

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_prepare_headmodel prepare_mesh_tetrahedral prepare_mesh_hexahedral

n = 71;

mri = [];
mri.dim = [n n n];
mri.tissue = zeros(mri.dim);
mri.transform = [
  1 0 0 -(n-1)/2-1
  0 1 0 -(n-1)/2-1
  0 0 1 -(n-1)/2-1
  0 0 0 1
  ];
mri.unit = 'mm';

[X, Y, Z] = ndgrid(1:n, 1:n, 1:n);
voxelpos = ft_warp_apply(mri.transform, [X(:) Y(:) Z(:)]);

ft_plot_mesh(voxelpos)
grid on
axis on

skin  = sqrt(sum(voxelpos.^2,2))<0.8*n/2;
skull = sqrt(sum(voxelpos.^2,2))<0.7*n/2;
brain = sqrt(sum(voxelpos.^2,2))<0.6*n/2;

mri.tissue(skin)  = 1;
mri.tissue(skull) = 2;
mri.tissue(brain) = 3;
mri.tissuelabel = {'skin', 'skull', 'brain'};

% mri.tissue(brain) = 1;
% mri.tissuelabel = {'brain'};

cfg = [];
cfg.funparameter = 'tissue';
ft_sourceplot(cfg, mri);

%%

cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = {'skin', 'skull', 'brain'};
cfg.numvertices = [100 200 300];
mesh1 = ft_prepare_mesh(cfg, mri);

figure
ft_plot_mesh(mesh1, 'facecolor', 'none', 'edgecolor', 'k')

%%
cfg = [];
cfg.method = 'iso2mesh';
cfg.tissue = {'skin', 'skull', 'brain'};
cfg.numvertices = [100 200 300];
mesh2 = ft_prepare_mesh(cfg, mri);

figure
ft_plot_mesh(mesh2, 'facecolor', 'none')

%%
cfg = [];
cfg.method = 'isosurface';
cfg.tissue = {'skin', 'skull', 'brain'};
cfg.numvertices = [100 200 300];
mesh3 = ft_prepare_mesh(cfg, mri);

figure
ft_plot_mesh(mesh3, 'facecolor', 'none', 'edgecolor', 'k')

%%
cfg = [];
cfg.method = 'hexahedral';
cfg.tissue = {'skin', 'skull', 'brain'};
mesh4 = ft_prepare_mesh(cfg, mri);

figure
ft_plot_mesh(mesh4, 'facecolor', 'none', 'edgecolor', 'k', 'surfaceonly', true)

%%
cfg = [];
cfg.method = 'tetrahedral';
cfg.tissue = {'skin', 'skull', 'brain'};
mesh5 = ft_prepare_mesh(cfg, mri);

figure
ft_plot_mesh(mesh5, 'facecolor', 'none', 'edgecolor', 'k', 'surfaceonly', true)

