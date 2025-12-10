function test_bug2959

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_sourceanalysis
% DATA private

load(dccnpath('/project/3031000.02/test/bug2959.mat')); 

source = ft_sourceanalysis(cfg, freq);
