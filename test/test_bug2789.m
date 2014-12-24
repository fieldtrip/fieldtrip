function test_bug2789

% WALLTIME 00:10:00
% MEM 1500mb

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2789.mat'));

% do not repeat the debugging here
cfg.debug = 'no';

% the following fails
ft_sourceanalysis(cfg, data);

