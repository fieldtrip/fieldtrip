function test_bug1527

% TEST test_bug1527
% TEST ft_sourceplot

cd home/common/matlab/fieldtrip/data/test/
load bug1527.mat

% this should be enough to reproduce the error
ft_sourceplot(cfg, source)
