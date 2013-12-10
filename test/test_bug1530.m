function test_bug1530

% MEM 2gb
% WALLTIME 00:03:08

% TEST test_bug1530
% TEST ft_sourceplot

% The problem: apparently ft_sourceplot fails on parameterselection
load('/home/common/matlab/fieldtrip/data/test/bug1530/cfg_sourceDiffIntNorm');
load('/home/common/matlab/fieldtrip/data/test/bug1530/sourceDiffIntNorm');

ft_sourceplot(cfg, sourceDiffIntNorm);
