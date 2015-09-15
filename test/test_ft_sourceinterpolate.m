function test_ft_sourceinterpolate

% MEM 4500mb
% WALLTIME 00:10:00

% TEST test_ft_sourceinterpolate
% TEST ft_sourceinterpolate

% See also test_bug2769 which goes over more interpolation options using fake data

% ft_sourceinterpolate is meant to interpolate functional data to anatomical data
% ft_sourceinterpolate(cfg, fun data, anat data)

% set plotting configuration used for all sections below:
cfgp = [];
cfgp.funparameter = 'pow';
cfgp.method       = 'slice';

%% 3D -> 3D
%  3D source to MRI 
mri = ft_read_mri('/home/common/matlab/fieldtrip/data/Subject01.mri');
%  get 'source'
% load('/home/common/matlab/fieldtrip/data/test/latest/source/meg/source_grid_timelock_trl_LCMV_keepnothing_ctf275.mat');
load('source3D')

cfg = [];
cfg.parameter = 'avg.pow';
i1 = ft_sourceinterpolate(cfg, source3D, mri);
ft_sourceplot(cfgp,i1);

%  3D source to 3D grid
% create 3D grid
[x,y,z]      = ndgrid(4:4:256,4:4:256,4:4:256);
grid.pos     = ft_warp_apply(mri.transform, [x(:) y(:) z(:)]); clear x y z
grid.dim     = [64 64 64];
grid.inside  = 1:size(grid.pos,1);
grid.outside = [];
grid.unit    = 'mm';

cfg = [];
cfg.parameter = 'avg.pow';
i2 = ft_sourceinterpolate(cfg, fun1, grid);
ft_sourceplot(cfgp,i2);


%% 2D > 2D
%  low resolution data -> high resolution grid

%  low resolution (single subject) data
% load('/project/3011020.09/MEG/A2020/coh/A2020_coh_sourcedata_delta_surface_ampnorm_sentpeak.mat');
load('fun2D');
fun2D.coh    = fun2D.avg.coh; % move to first level

%   high resolution (single subject) mesh
% load('/project/3011020.09/MEG/A2020/anatomy/A2020sourcemodel2Dsurfreg');
load('surf2D_hires');

% Interpolation with 'nearest'
cfg   = [];
cfg.parameter = 'coh';
cfg.interpmethod = 'nearest';   % default for 2D mesh
i3a = ft_sourceinterpolate(cfg, fun2D, surf2D_hires.orig);
cfgp.method       = 'surface';
cfgp.funparameter = 'coh';
ft_sourceplot(cfgp,i3a);


% Interpolation with 'smudge'
% FIXME:  error 'the smudge method needs triangle definition')
cfg.interpmethod = 'smudge';   % option uses sphereradius
i3b = ft_sourceinterpolate(cfg, fun2D, surf2D_hires.orig);
ft_sourceplot(cfgp,i3b)

%% 2D > 3D
% 2D surface to MRI
% load('/home/common/matlab/fieldtrip/data/test/latest/source/meg/source_sheet_timelock_trl_MNE_keepnothing_ctf275.mat');
load('fun2D_MNE')

cfg = [];
cfg.parameter = 'avg.pow';
cfg.downsample = 4; % serious downsampling is needed to keep it fitting in memory with 600 time points
i4a = ft_sourceinterpolate(cfg, fun2D_MNE, mri);

% FIXME: error: no funparameter no anaparameter
cfgp.method       = 'slice';
cfgp.funparameter = 'pow';
ft_sourceplot(cfgp,i4a)

% 2D source to 3D grid
cfg = [];
cfg.parameter = 'avg.pow';
cfg.downsample = 4; % serious downsampling is needed to keep it fitting in memory with 600 time points
i4b = ft_sourceinterpolate(cfg, source2, grid);


%% 3D > 2D
 %  volume data to cortical mesh 
 
 
