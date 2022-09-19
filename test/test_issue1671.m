function test_issue1671

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_megplanar ft_apply_montage

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/connectivity/data.mat'));

cfg = [];
cfg.method = 'template';
cfg.template = 'ctf151_neighb.mat';
neighb = ft_prepare_neighbours(cfg);

cfg = [];
cfg.method = 'sincos';
cfg.neighbours = neighb;
datplan = ft_megplanar(cfg, data);

assert(any(strcmp(data.label, 'EMGlft')));
assert(any(strcmp(datplan.label, 'EMGlft')));

end