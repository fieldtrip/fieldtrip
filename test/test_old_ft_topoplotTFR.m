function test_old_ft_topoplotTFR

% MEM 1gb
% WALLTIME 00:20:00


% TEST_FT_TOPOPLOTTFR
% This script tests the ft_topoplotTFR function and should display a figure
% with the ctf275 layout showing power decreases at the parietal lobes

% load time-frequency data
datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/');
fprintf('loading data\n');
%load observe_comm_moves_freqmtmconvol.mat
load(fullfile(datadir,'observe_comm_moves_freqmtmconvol.mat'));

% topoplot a low frequency band
figure;
cfg              = [];
cfg.baselinetype = 'relchange';
cfg.baseline     = [-2 -.5];
cfg.xlim         = [.5 3];
cfg.ylim         = [8 12];
cfg.zlim         = [-1 1];
ft_topoplotTFR(cfg, obs_lo);

% topoplot a high frequency band
figure;
cfg              = [];
cfg.baselinetype = 'relchange';
cfg.baseline     = [-2 -.5];
cfg.xlim         = [.5 3];
cfg.ylim         = [30 40];
cfg.zlim         = [-1 1];
ft_topoplotTFR(cfg, obs_hi);
