function test_bug1212

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_layoutplot

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1212.mat'));

cfg = [];
cfg.grad = hdr.grad;
ft_layoutplot(cfg, [])
