function test_bug2399

% MEM 1gb
% WALLTIME 00:20:00
% DEPENDENCY ft_sourceanalysis ft_prepare_vol_sens channelposition
% DATA private

load(dccnpath('/project/3031000.02/test/bug2399.mat'));
 
ft_sourceanalysis(cfg, tl);
