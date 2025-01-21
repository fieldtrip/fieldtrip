function test_issue2193

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_apply_montage ft_prepare_montage 
% DATA private

load(dccnpath('/project/3031000.02/test/issue2193.mat'));
data = ft_preprocessing(cfg, data_segmented);
