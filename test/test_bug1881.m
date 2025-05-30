function test_bug1881

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_selectdata
% DATA private

filename = dccnpath('/project/3031000.02/test/bug1881.mat');
load(filename);

cfg        = [];
cfg.foilim = [0 90];
output     = ft_selectdata(cfg, freq);

assert(length(output.freq)==419, 'incorrect frequency selection')
