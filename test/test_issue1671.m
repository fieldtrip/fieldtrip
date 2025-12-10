function test_issue1671

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_megplanar ft_apply_montage
% DATA public

load(dccnpath('/project/3031000.02/external/download/tutorial/connectivity/data.mat'));

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