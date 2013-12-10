function test_bug2399

% WALLTIME 0:05:00
% MEM 2000mb

% TEST test_bug2399
% TEST ft_sourceanalysis ft_prepare_vol_sens channelposition

load(dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2399.mat'));
 
ft_sourceanalysis(cfg, tl);
