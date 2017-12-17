function test_bug1051

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_neighbours ft_megplanar ft_apply_montage

% the bug consists of ft_apply_montage not adequately dealing with 
% sensor descriptions that contain coilori/pos chanori/pos.
% The chanori/chanpos get lost along the way

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg'));
load preproc_ctf151

cfg.method = 'triangulation';
neighbours = ft_prepare_neighbours(cfg, data);
cfg        = [];
cfg.planarmethod = 'sincos';
cfg.neighbours   = neighbours;
planar     = ft_megplanar(cfg, data);

% the loss of chanpos seems to be appropriate, but the subsequent call
% to e.g. freqanalysis fails when fixsens applies a too strict checking
% relaxing the conditional statements prevents the crash

cfg = [];
cfg.method = 'mtmfft';
cfg.foilim = [0 20];
cfg.taper = 'hanning';
freq = ft_freqanalysis(cfg, planar);
