function test_bug1881

% WALLTIME 00:10:00
% MEM 1024mb

% TEST ft_selectdata

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1881.mat');
load(filename);

cfg        = [];
cfg.foilim = [0 90];
output     = ft_selectdata(cfg, freq);

assert(length(output.freq)==419, 'incorrect frequency selection')
