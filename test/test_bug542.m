function test_bug542

% TEST test_bug542
% TEST ft_multiplotER

cd(dccnfilename('/home/common/matlab/fieldtrip/data/test'));
load bug542.mat

stopwatch = tic;
ft_multiplotER(cfg, att_dep_ipsi);
toc(stopwatch)

