function test_ft_defacevolume

% MEM 8gb
% WALLTIME 00:10:00
% DEPENDENCY ft_defacevolume

mri = [];
mri.anatomy = randn(181,217,181);
mri.dim = [181 217 181];
mri.transform = eye(4);
mri.unit = 'cm';
mri.coordsys = 'ctf';

cfg = [];
cfg.method = 'spm';
mriout = ft_defacevolume(cfg, mri);
