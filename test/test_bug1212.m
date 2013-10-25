function test_bug1212

% MEM 1gb
% WALLTIME 00:03:04

% TEST test_bug1212
% TEST ft_layoutplot

load /home/common/matlab/fieldtrip/data/test/bug1212.mat

cfg = [];
cfg.grad = hdr.grad;
ft_layoutplot(cfg, [])
