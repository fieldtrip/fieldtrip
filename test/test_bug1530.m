function test_bug1530

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_sourceplot

% The problem: apparently ft_sourceplot fails on parameterselection
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1530/cfg_sourceDiffIntNorm'));
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1530/sourceDiffIntNorm'));

% update surffile according to new filename
cfg.surffile = 'surface_inflated_both.mat';
ft_sourceplot(cfg, sourceDiffIntNorm);
