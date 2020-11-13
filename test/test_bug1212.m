function test_bug1212

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_layoutplot

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1212.mat'));

cfg = [];
cfg.grad = hdr.grad;
ft_layoutplot(cfg, [])
