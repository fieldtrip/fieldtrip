function test_pull1427

% MEM 6gb
% WALLTIME 1:00:00
% DEPENDENCY ft_prepare_mesh ft_prepare_headmodel

%%

% This function creates an hexahedral and tetrahedral volumetric mesh from a
% three-compartments volume to be passed as input to ft_prepare_headmodel with
% the method 'simbio'.

segprob = [];
segprob.brain = false(10,10,10); segprob.brain(4:7,4:7,4:7) = true;
segprob.skull = false(10,10,10); segprob.skull(3:8,3:8,3:8) = true;
segprob.scalp = false(10,10,10); segprob.scalp(2:9,2:9,2:9) = true;
% segprob.air = true(10,10,10);
% segprob.air(2:9,2:9,2:9) = false;
%segprob.brain = false(11,11,11); segprob.brain(4:8,4:8,4:8) = true;
%segprob.skull = false(11,11,11); segprob.skull(3:9,3:9,3:9) = true;
%segprob.scalp = false(11,11,11); segprob.scalp(2:10,2:10,2:10) = true;
segprob.dim = size(segprob.brain);
segprob.unit = 'cm';
segprob.coordsys = 'ctf';
segprob.transform = eye(4);
segprob.transform(1,4) = -0.5;
segprob.transform(2,4) = -0.5;
segprob.transform(3,4) = -0.5;

% it is more difficult to visualize a probabilistic segmentation than an indexed one
segindx = ft_datatype_segmentation(segprob, 'segmentationstyle', 'indexed');

cfg = [];
cfg.funparameter = 'tissue';
cfg.method = 'ortho';
cfg.location = [5 5 5]; % this is the center of the volume, in this plot it will be rounded off to the nearest voxel
ft_sourceplot(cfg, segindx);

% determine the range of the bounding box
[X, Y, Z] = ndgrid(1:segprob.dim(1), 1:segprob.dim(2), 1:segprob.dim(3));
voxpos = ft_warp_apply(segprob.transform, [X(:) Y(:) Z(:)]);
minmaxpos(:,1) = min(voxpos) - 0.5;
minmaxpos(:,2) = max(voxpos) + 0.5;

%%

cfg = [];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
mesh_vol_hex1 = ft_prepare_mesh(cfg, segprob);

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
hold on
ft_plot_mesh(mesh_vol_hex1, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [0.33, 0.43, 0.53]; % order follows mesh.tissuelabel
ft_prepare_headmodel(cfg, mesh_vol_hex1)

%%

% make a volume with 2x the resolution
cfg = [];
cfg.resolution = 1/2;
cfg.xrange     = minmaxpos(1,:); % this is more robust than the dim, in case the origin is not inside the volume
cfg.yrange     = minmaxpos(2,:);
cfg.zrange     = minmaxpos(3,:);
cfg.method     = 'nearest';
segprob2 = ft_volumereslice(cfg, segprob);
segindx2 = ft_volumereslice(cfg, segindx);

cfg = [];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
mesh_vol_hex2 = ft_prepare_mesh(cfg, segprob2);

% this one looks OK in terms of alignment
figure
% ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
ft_plot_ortho(segindx2.tissue, 'transform', segindx2.transform, 'location', [5 5 5], 'style', 'intersect'); % this one looks OK in terms of alignment
hold on
ft_plot_mesh(mesh_vol_hex2, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [0.33, 0.43, 0.53]; % order follows mesh.tissuelabel
ft_prepare_headmodel(cfg, mesh_vol_hex1)

%%

% make a volume with 3x the resolution
cfg = [];
cfg.resolution = 1/3;
cfg.xrange     = minmaxpos(1,:); % this is more robust than the dim, in case the origin is not inside the volume
cfg.yrange     = minmaxpos(2,:);
cfg.zrange     = minmaxpos(3,:);
cfg.method     = 'nearest';
segprob3 = ft_volumereslice(cfg, segprob);
segindx3 = ft_volumereslice(cfg, segindx);

cfg = [];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
mesh_vol_hex3 = ft_prepare_mesh(cfg, segprob3);

figure
% ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
% ft_plot_ortho(segindx2.tissue, 'transform', segindx2.transform, 'location', [5 5 5], 'style', 'intersect');
ft_plot_ortho(segindx3.tissue, 'transform', segindx3.transform, 'location', [5 5 5], 'style', 'intersect'); % this one looks OK in terms of alignment
hold on
ft_plot_mesh(mesh_vol_hex3, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [0.33, 0.43, 0.53]; % order follows mesh.tissuelabel
ft_prepare_headmodel(cfg, mesh_vol_hex3)

%%

cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segprob);

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
hold on
ft_plot_mesh(mesh_vol_tet, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [0.33, 0.43, 0.53]; % order follows mesh.tissuelabel
ft_prepare_headmodel(cfg, mesh_vol_tet)
