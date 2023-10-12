function test_bug1035

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_mulitplotER ft_prepare_layout
% DATA private

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1035.mat

ft_multiplotER(cfg, avg151)
