function test_bug2544


% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_sourcegrandaverage getdimord ft_selectdata
% DATA private

load(dccnpath('/project/3031000.02/test/bug2544.mat'));

cfg = [];
grandavg = ft_sourcegrandaverage(cfg, tmp{:});

