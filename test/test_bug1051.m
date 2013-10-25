function test_bug1051

% MEM 1gb
% WALLTIME 00:03:03

% TEST test_bug1051 ft_megplanar ft_apply_montage

% the bug consists of ft_apply_montage not adequately dealing with 
% sensor descriptions that contain coilori/pos chanori/pos.
% The chanori/chanpos get lost along the way

cd /home/common/matlab/fieldtrip/data/test/latest/raw/meg/
load preproc_ctf151

cfg.method = 'triangulation';
neighbours = ft_prepare_neighbours(cfg, data);
cfg        = [];
cfg.planarmethod = 'sincos';
cfg.neighbours   = neighbours;
planar     = ft_megplanar(cfg, data);
