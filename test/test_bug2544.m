function test_bug2544


% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_sourcegrandaverage getdimord ft_selectdata

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2544.mat'));

cfg = [];
grandavg = ft_sourcegrandaverage(cfg, tmp{:});

