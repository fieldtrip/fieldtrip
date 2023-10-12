function test_bug1212

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_layoutplot
% DATA private

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1212.mat'));

cfg = [];
cfg.grad = hdr.grad;
ft_layoutplot(cfg, [])
