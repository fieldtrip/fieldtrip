function failed_bug2769

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_bug2769
% TEST ft_sourceplot ft_sourceinterpolate

clear all
close all

%% construct initial volume and surface

mri1 = [];
mri1.dim = [11 11 11];
mri1.transform = [
  1 0 0 -6
  0 1 0 -6
  0 0 1 -6
  0 0 0  1
  ];

[X, Y, Z] = ndgrid(1:mri1.dim(1), 1:mri1.dim(2), 1:mri1.dim(3));
pos1 = ft_warp_apply(mri1.transform, [X(:) Y(:) Z(:)]);
seg1 = sqrt(sum(pos1.^2, 2))<4;

mri1.anatomy = reshape(seg1, mri1.dim);
mri1.seg     = reshape(seg1, mri1.dim);


mri2 = [];
mri2.dim = [19 19 19];
mri2.transform = [
  1 0 0 -10
  0 1 0 -10
  0 0 1 -10
  0 0 0  1
  ];

[X, Y, Z] = ndgrid(1:mri2.dim(1), 1:mri2.dim(2), 1:mri2.dim(3));
pos2 = ft_warp_apply(mri2.transform, [X(:) Y(:) Z(:)]);
seg2 = sqrt(sum(pos2.^2, 2))<4;

mri2.anatomy = reshape(seg2, mri2.dim);
mri2.seg     = reshape(seg2, mri2.dim);

% add some noise, otherwise it resembles a segmentation
mri1.anatomy = mri1.anatomy + 0.1*randn(size(mri1.anatomy));
mri2.anatomy = mri2.anatomy + 0.1*randn(size(mri2.anatomy));

cfg = [];
ft_sourceplot(cfg, mri1);
ft_sourceplot(cfg, mri2);

volume1 = mri1;
volume1.pow = reshape(pos1(:,3), mri1.dim);

volume2 = mri2;
volume2.pow = reshape(pos2(:,3), mri2.dim);

cfg = [];
cfg.tissue = 'seg';
cfg.numvertices = 1000;
mesh1 = ft_prepare_mesh(cfg, mri1);

cfg = [];
cfg.tissue = 'seg';
cfg.numvertices = 1500;
mesh2 = ft_prepare_mesh(cfg, mri2);

surface1 = mesh1;
surface1.pow = surface1.pnt(:,3);

surface2 = mesh2;
surface2.pow = surface2.pnt(:,3);

%% interpolate volume data onto surface

cfg = [];
cfg.parameter = 'pow';
interpA = ft_sourceinterpolate(cfg, volume1, surface1);

assert(ft_datatype(interpA, 'source'));


%% interpolate surface data onto volume

cfg = [];
cfg.parameter = 'pow';
interpB = ft_sourceinterpolate(cfg, surface1, volume1);

assert(ft_datatype(interpB, 'volume'));

%% interpolate volume data onto another volume

cfg = [];
cfg.parameter = 'pow';
interpC = ft_sourceinterpolate(cfg, volume1, mri2);

assert(ft_datatype(interpC, 'volume'));

%% interpolate surface data onto another surface

cfg = [];
cfg.parameter = 'pow';
interpD = ft_sourceinterpolate(cfg, surface1, mesh2);

assert(ft_datatype(interpD, 'source'));

%%
cfg = [];
cfg.funparameter = 'pow';
cfg.method = 'surface';

ft_sourceplot(cfg, interpA);
ft_sourceplot(cfg, volume1, surface1);

cfg.method = 'ortho';
ft_sourceplot(cfg, interpB);
ft_sourceplot(cfg, surface1, volume1);

cfg.method = 'ortho';
ft_sourceplot(cfg, interpC);
ft_sourceplot(cfg, volume1, mri2);

cfg.method = 'surface';
ft_sourceplot(cfg, interpD);
ft_sourceplot(cfg, surface1, mesh2);



