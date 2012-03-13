function test_bug1212

% TEST test_bug1212
% TEST ft_layoutplot

load /home/common/matlab/fieldtrip/data/test/bug1212.mat

cfg = [];
cfg.grad = hdr.grad;
ft_layoutplot(cfg, [])
