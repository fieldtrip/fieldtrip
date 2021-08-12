function test_ft_volumebiascorrect

% MEM 8gb
% WALLTIME 00:60:00
% DEPENDENCY ft_volumebiascorrect

mri = [];
mri.anatomy = randn(181,217,181);
mri.dim = [181 217 181];
mri.transform = eye(4);
mri.unit = 'cm';
mri.coordsys = 'ctf';

cfg = [];
mriout = ft_volumebiascorrect(cfg, mri);
