function test_bug1043

% MEM 1gb
% WALLTIME 00:03:08

% TEST test_bug1043 
% TEST ft_megplanar ft_apply_montage yokogawa2grad channelposition

% the bug consists of ft_apply_montage not adequately dealing with
% sensor descriptions that contain coilori/pos chanori/pos.
% The chanori/chanpos get lost along the way

cd /home/common/matlab/fieldtrip/data/test/latest/raw/meg/
load preproc_ctf151

cfg        = [];
cfg.method = 'triangulation';
neighbours = ft_prepare_neighbours(cfg, data);

cfg               = [];
cfg.planarmethod  = 'sincos';
cfg.neighbours    = neighbours;
planar            = ft_megplanar(cfg, data);

% this is another test case that failed to work according to the bug report

cd /home/common/matlab/fieldtrip/data/test/bug1043
cfg = [];
cfg.dataset = 'chki250110.0100.fs1khz.corr_ch_names-ave.ave';
data = ft_preprocessing(cfg);
