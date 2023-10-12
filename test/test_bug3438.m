function test_bug3438

% MEM 6gb
% WALLTIME 00:20:00
% DEPENDENCY ft_read_headshape
% DATA private

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3438'))

% the following failed prior to fixing this bug
shape = ft_read_headshape('OBJ.obj');
shape = ft_read_headshape('PLY.ply');
