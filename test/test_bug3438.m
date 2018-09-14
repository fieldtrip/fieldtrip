function test_bug3438

% WALLTIME 00:20:00
% MEM 5gb

% TEST ft_read_headshape

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3438'))

% the following failed prior to fixing this bug
shape = ft_read_headshape('OBJ.obj');
shape = ft_read_headshape('PLY.ply');
