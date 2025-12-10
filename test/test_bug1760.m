function test_bug1760

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_multiplotER ft_multiplotTFR
% DATA private

cd(dccnpath('/project/3031000.02/test'))
load bug1760.mat

figure
ft_multiplotER(cfg, freq);

