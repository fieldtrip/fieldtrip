function test_bug542

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_multiplotER

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug542.mat'));

% this is confusing the plotting function
att_dep_ipsi = rmfield(att_dep_ipsi, 'time');

stopwatch = tic;
ft_multiplotER(cfg, att_dep_ipsi);
toc(stopwatch)
