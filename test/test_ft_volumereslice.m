function test_ft_volumereslice

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_volumereslice SPM

mri = [];
mri.anatomy = randn(181,217,181);
mri.dim = [181 217 181];
mri.transform = eye(4);
mri.unit = 'cm';
mri.coordsys = 'ctf';

cfg = [];
cfg.method = 'linear';
cfg.resolution = 2;
cfg.dim = [181 217 181];
mriout = ft_volumereslice(cfg, mri);