function inspect_bug3013

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_sourceplot ft_plot_ortho
% DATA no

%%

mri = [];
mri.anatomy = randn(10,20,40);
mri.dim = size(mri.anatomy);
mri.transform = [
  1 0 0 0
  0 1 0 0
  0 0 1 0
  0 0 0 1
  ];

%%

cfg = [];
ft_sourceplot(cfg, mri); % use defaults
drawnow

cfg.voxelratio = 'square';

cfg.axisratio = 'square';
ft_sourceplot(cfg, mri);
drawnow
cfg.axisratio = 'voxel';
ft_sourceplot(cfg, mri);
drawnow
cfg.axisratio = 'data';
ft_sourceplot(cfg, mri);
drawnow

cfg.voxelratio = 'data';

cfg.axisratio = 'square';
ft_sourceplot(cfg, mri);
drawnow
cfg.axisratio = 'voxel';
ft_sourceplot(cfg, mri);
drawnow
cfg.axisratio = 'data';
ft_sourceplot(cfg, mri);
drawnow


%%

mri = [];
mri.anatomy = randn(10,20,40);
mri.dim = size(mri.anatomy);
mri.transform = [
  1/10 0 0 0
  0 1/20 0 0
  0 0 1/40 0
  0 0 0 1
  ];

%%

cfg = [];
ft_sourceplot(cfg, mri); % use defaults
drawnow

cfg.voxelratio = 'square';

cfg.axisratio = 'square';
ft_sourceplot(cfg, mri);
drawnow
cfg.axisratio = 'voxel';
ft_sourceplot(cfg, mri);
drawnow
cfg.axisratio = 'data';
ft_sourceplot(cfg, mri);
drawnow

cfg.voxelratio = 'data';

cfg.axisratio = 'square';
ft_sourceplot(cfg, mri);
drawnow
cfg.axisratio = 'voxel';
ft_sourceplot(cfg, mri);
drawnow
cfg.axisratio = 'data';
ft_sourceplot(cfg, mri);
drawnow
