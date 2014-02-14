function test_bug2464

% MEM 150mb
% WALLTIME 00:10:00

% TEST test_bug2464
% TEST ft_selectdata

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2464.mat');
load(filename);

cfg = [];
cfg.foilim = [10 20];
sel = ft_selectdata(cfg, freq);