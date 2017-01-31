function failed_bug2850

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_sourceanalysis ft_eloreta

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

load(dccnpath('/home/common/matlab/fieldtrip/data/test/avgFIC.mat'));

vol = [];
vol.r = 10;
vol.o = [0 0 4];
vol.unit = 'cm';

cfg = [];
cfg.vol = vol;
cfg.grad = avgFIC.grad;
cfg.grid.resolution = 1;
cfg.grid.unit = 'cm';
grid = ft_prepare_sourcemodel(cfg);

% due to the spherical volume conductor and the regular grid, some dipole
% positions would get a leadfield of rank 1, which eloreta does not like
grid.pos = grid.pos + 0.01*randn(size(grid.pos));

cfg = [];
cfg.vol = vol;
cfg.grid = grid;
cfg.channel = 'MEG';
grid = ft_prepare_leadfield(cfg, avgFIC);

cfg = [];
cfg.vol = vol;
cfg.grid = grid;
cfg.method = 'eloreta';
cfg.eloreta.keepleadfield = 'yes';
source = ft_sourceanalysis(cfg, avgFIC);

assert(isfield(source, 'avg'))
assert(isfield(source.avg, 'mom'))
assert(isfield(source.avg, 'pow'))
assert(~isfield(source.avg, 'inside'))    % this should be at the top level
assert(~isfield(source.avg, 'leadfield')) % this should be at the top level

insideindx = find(source.inside);
assert(size(source.avg.mom{insideindx(1)},1)==3);       % this should be 3*Ntime
assert(size(source.avg.leadfield{insideindx(1)},2)==3); % this should be Nchan*3


