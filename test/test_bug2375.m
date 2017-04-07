function test_bug2375

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_headmodel ft_headmodel_localspheres

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2375/localspheres_bug.mat'));

vol = ft_prepare_headmodel(cfg, headshape);

% they should all have the same sphere, since the cfg is faulty
assert(all(vol.r==vol.r(1)));

% fix the cfg and try again
cfg = rmfield(cfg, 'radius');
cfg = rmfield(cfg, 'maxradius');
vol = ft_prepare_headmodel(cfg, headshape);

% they should not all have the same sphere any more
assert(~all(vol.r==vol.r(1)));
