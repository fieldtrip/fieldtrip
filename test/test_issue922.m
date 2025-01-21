function test_issue922

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqanalysis
% DATA private

%test case for ft_freqanalysis rounding bug

load(dccnpath('/project/3031000.02/test/issue922.mat'));
[freq] = ft_freqanalysis(cfg, test);
