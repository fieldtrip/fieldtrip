function test_ft_sourceinterpolate

% MEM 4500mb
% WALLTIME 00:10:00

% TEST test_ft_sourceinterpolate
% TEST ft_sourceinterpolate

% See also test_bug2769 which goes over more interpolation options, but uses fake data

% 3D source to MRI
mri = ft_read_mri('/home/common/matlab/fieldtrip/data/Subject01.mri');
load('/home/common/matlab/fieldtrip/data/test/latest/source/meg/source_grid_timelock_trl_LCMV_keepnothing_ctf275.mat');
source1 = source;

cfg = [];
cfg.parameter = 'avg.pow';
i1 = ft_sourceinterpolate(cfg, source1, mri);

% 3D source to 3D grid
[x,y,z] = ndgrid(4:4:256,4:4:256,4:4:256);
grid.pos = ft_warp_apply(mri.transform, [x(:) y(:) z(:)]); clear x y z
grid.dim = [64 64 64];
grid.inside = 1:size(grid.pos,1);
grid.outside = [];
grid.unit = 'mm';

cfg = [];
cfg.parameter = 'avg.pow';
i2 = ft_sourceinterpolate(cfg, source1, grid);

% 2D source to MRI
load('/home/common/matlab/fieldtrip/data/test/latest/source/meg/source_sheet_timelock_trl_MNE_keepnothing_ctf275.mat');
source2 = source;
source2.unit = 'cm'; % FIXME ft_convert_units determines source2 to have 'dm' as unit?

cfg = [];
cfg.parameter = 'avg.pow';
cfg.downsample = 4; % serious downsampling is needed to keep it fitting in memory with 600 time points
i3 = ft_sourceinterpolate(cfg, source2, mri);

% 2D source to 3D grid
cfg = [];
cfg.parameter = 'avg.pow';
cfg.downsample = 4; % serious downsampling is needed to keep it fitting in memory with 600 time points
i4 = ft_sourceinterpolate(cfg, source2, grid);
