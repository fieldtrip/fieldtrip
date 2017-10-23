function test_bug2399

% MEM 2000mb
% WALLTIME 00:20:00

% TEST ft_sourceanalysis ft_prepare_vol_sens channelposition

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2399.mat'));
 
ft_sourceanalysis(cfg, tl);
