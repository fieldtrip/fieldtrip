function test_issue2193

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_apply_montage ft_prepare_montage 

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue2193.mat'));
data = ft_preprocessing(cfg, data_segmented);
