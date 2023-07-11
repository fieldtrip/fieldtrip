function test_issue2274

% WALLTIME 00:10:00
% MEM 4gb
% DEPENDENCY GETDIMORD FT_SOURCEPLOT

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue2274.mat'));

ft_sourceplot(cfg, sourceProj)
