function test_bug2959

% WALLTIME 00:10:00
% MEM 2000mb

% TEST ft_sourceanalysis

load(fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test'),'bug2959.mat')); 

source = ft_sourceanalysis(cfg, freq);
