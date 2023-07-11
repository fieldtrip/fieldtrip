function failed_bug2850

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_sourceanalysis ft_inverse_eloreta

load(dccnpath('/home/common/matlab/fieldtrip/data/test/avgFIC.mat'));

vol = [];
vol.r = 10;
vol.o = [0 0 4];
vol.unit = 'cm';

cfg = [];
cfg.headmodel = vol;
cfg.grad = avgFIC.grad;
cfg.resolution = 1;
cfg.unit = 'cm';
grid = ft_prepare_sourcemodel(cfg);

% due to the spherical volume conductor and the regular grid, some dipole
% positions would get a leadfield of rank 1, which eloreta does not like
sourcemodel.pos = sourcemodel.pos + 0.01*randn(size(sourcemodel.pos));

cfg = [];
cfg.headmodel = vol;
cfg.sourcemodel = grid;
cfg.channel = 'MEG';
grid = ft_prepare_leadfield(cfg, avgFIC);

cfg = [];
cfg.headmodel = vol;
cfg.sourcemodel = grid;
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


