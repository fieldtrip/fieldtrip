function test_bug1530

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_sourceplot
% DATA private

% The problem: apparently ft_sourceplot fails on parameterselection
load(dccnpath('/project/3031000.02/test/bug1530/cfg_sourceDiffIntNorm.mat'));
load(dccnpath('/project/3031000.02/test/bug1530/sourceDiffIntNorm.mat'));

% update surffile according to new filename
cfg.surffile = 'surface_inflated_both.mat';
ft_sourceplot(cfg, sourceDiffIntNorm);
