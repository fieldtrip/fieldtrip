function test_bug542

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_multiplotER

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'));
load bug542.mat

stopwatch = tic;
ft_multiplotER(cfg, att_dep_ipsi);
toc(stopwatch)

