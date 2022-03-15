function test_issue922

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqanalysis

%test case for ft_freqanalysis rounding bug

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue922.mat'));
[freq] = ft_freqanalysis(cfg, test);
