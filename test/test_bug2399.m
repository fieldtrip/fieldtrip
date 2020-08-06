function test_bug2399

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_sourceanalysis ft_prepare_vol_sens channelposition

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2399.mat'));
 
ft_sourceanalysis(cfg, tl);
