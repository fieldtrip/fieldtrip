function test_ft_sourceinterpolate

% WALLTIME 00:03:14

% TEST test_ft_sourceinterpolate
% TEST ft_sourceinterpolate

% 3D source to MRI
mri = ft_read_mri('/home/common/matlab/fieldtrip/data/Subject01.mri');
load('/home/common/matlab/fieldtrip/data/test/latest/source/meg/source_grid_timelock_trl_LCMV_keepnothing_ctf275.mat');
source1 = source;

cfg = [];
cfg.parameter = 'avg.pow';
i1 = ft_sourceinterpolate(cfg, source1, mri);

% 3D source to 3D grid
[x,y,z] = ndgrid(4:4:256,4:4:256,4:4:256);
grid.pos = warp_apply(mri.transform, [x(:) y(:) z(:)]); clear x y z
grid.dim = [64 64 64];
grid.inside = 1:size(grid.pos,1);
grid.outside = [];
grid.unit = 'mm';

i2 = ft_sourceinterpolate(cfg, source1, grid);

% 3D source to 2D sheet

% 2D source to MRI
load('/home/common/matlab/fieldtrip/data/test/latest/source/meg/source_sheet_timelock_trl_MNE_keepnothing_ctf275.mat');
source2 = source;

% FIXME ft_convert_units determines source2 to have 'dm' as unit?
source2.unit = 'cm';

% i4 = ft_sourceinterpolate(cfg, source2, mri);

% 2D source to 3D grid
% i5 = ft_sourceinterpolate(cfg, source2, grid);
