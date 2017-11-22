function test_bug1035

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_mulitplotER ft_prepare_layout

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1035.mat

ft_multiplotER(cfg, avg151)


