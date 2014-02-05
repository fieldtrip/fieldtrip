function test_bug1894

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1894
% TEST ft_singleplotTFR ft_daattype_freq ft_datatype_sens ft_chantype

if ispc
  load h:\common\matlab\fieldtrip\data\test\bug1894.mat
else
  load /home/common/matlab/fieldtrip/data/test/bug1894.mat
end

cfg = [];
ft_multiplotTFR (cfg, freq);
ft_singleplotTFR(cfg, freq);
ft_topoplotTFR  (cfg, freq);

ft_datatype_sens(freq.grad);
