function test_bug1035

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_mulitplotER ft_prepare_layout
% DATA private

cd(dccnpath('/project/3031000.02/test'))
load bug1035.mat

ft_multiplotER(cfg, avg151)
