function test_bug1288

% this function serves to create planar gradient data and combined planar
% gradient data, for testing purposes of a fix for bug 1288

load('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat');

cfg = [];
cfg.method = 'distance';
cfg.grad   = data.grad;
neighbours = ft_prepare_neighbours(cfg);

cfg = [];
cfg.planarmethod = 'sincos';
cfg.neighbours   = neighbours;
datap = ft_megplanar(cfg, data);

cfg = [];
datac = ft_combineplanar(cfg, datap);

cfg = [];
cfg.grad = datac.grad;
lay = ft_prepare_layout(cfg);

% seems to run through without error