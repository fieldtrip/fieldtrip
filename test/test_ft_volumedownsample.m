function test_ft_volumedownsample

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_volumedownsample SPM

mri = [];
mri.anatomy = randn(181,217,181);
mri.dim = [181 217 181];
mri.transform = eye(4);
mri.unit = 'cm';
mri.coordsys = 'ctf';

cfg = [];
cfg.downsample = 2;
cfg.parameter = 'anatomy';
mriout = ft_volumedownsample(cfg, mri);