function test_pull1427

% MEM 6gb
% WALLTIME 1:00:00
% DEPENDENCY ft_prepare_mesh ft_prepare_headmodel

%%

% This function creates an hexahedral and tetrahedral volumetric mesh from a
% two-compartments volume to be passed as input to ft_prepare_headmodel with
% the method 'simbio'.

segprob = [];
segprob.brain = false(10,10,10); segprob.brain(4:7,4:7,4:7) = true;
segprob.skull = false(10,10,10); segprob.skull(3:8,3:8,3:8) = true;
segprob.scalp = false(10,10,10); segprob.scalp(2:9,2:9,2:9) = true;
segprob.dim = size(segprob.brain);
segprob.unit = 'mm';
segprob.coordsys = 'ctf';
segprob.transform = eye(4);

% it is more difficult to visualize a probabilistic segmentation than an indexed one
segindx = ft_datatype_segmentation(segprob, 'segmentationstyle', 'indexed');

cfg = [];
cfg.funparameter = 'seg';
cfg.method = 'ortho';
cfg.location = [5.5 5.5 5.5]; % this is the center of the volume, in this plot it will be rounded off to the nearest voxel
ft_sourceplot(cfg, segindx);

%%

cfg = [];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
cfg.resolution = 1;
mesh_vol_hex = ft_prepare_mesh(cfg, segprob);

figure
ft_plot_ortho(segindx.seg, 'transform', segindx.transform, 'location', [5.5 5.5 5.5], 'style', 'intersect');
hold on
ft_plot_mesh(mesh_vol_hex, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [0.33, 0.43, 0.53]; % order follows mesh.tissuelabel
ft_prepare_headmodel(cfg, mesh_vol_hex)

%%

cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segprob);

figure
ft_plot_ortho(segindx.seg, 'transform', segindx.transform, 'location', [5.5 5.5 5.5], 'style', 'intersect');
hold on
ft_plot_mesh(mesh_vol_tet, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [0.33, 0.43, 0.53]; % order follows mesh.tissuelabel
ft_prepare_headmodel(cfg, mesh_vol_tet)
