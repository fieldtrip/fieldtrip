function test_bug1894

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_singleplotTFR ft_daattype_freq ft_datatype_sens ft_chantype

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1894.mat'));

cfg = [];
ft_multiplotTFR (cfg, freq);
ft_singleplotTFR(cfg, freq);
ft_topoplotTFR  (cfg, freq);

ft_datatype_sens(freq.grad);
