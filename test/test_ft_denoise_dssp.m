function test_ft_denoise_dssp

% MEM 3gb
% WALLTIME 00:10:00
% DEPENDENCY ft_denoise_dssp

% below is based on simulated data using Kensuke's code, and should behave as the
% bdssp_main script that he uses to demonstrate the algorithm
load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull800/data.mat'));
tmp = load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull800/grid.mat'));

cfg = [];
cfg.sourcemodel = tmp.grid;
cfg.dssp.n_space = 30;
cfg.dssp.n_in    = 30;
cfg.dssp.n_out   = 30;
dataout = ft_denoise_dssp(cfg, data);

% the lowsnr data, has multiple trials, in which the 'blip' of activity has
% a slightly different latency each trial.
load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull800/data_lowsnr.mat'));
dataout2 = ft_denoise_dssp(cfg, data);

cfg = [];
cfg.method = 'mtmconvol';
cfg.foi = 0:2:20;
cfg.toi = data.time{1}(1:15:end);
cfg.t_ftimwin = ones(1,numel(cfg.foi))./2;
cfg.tapsmofrq = ones(1,numel(cfg.foi)).*2;
freq = ft_freqanalysis(cfg, data);
freq2 = ft_freqanalysis(cfg, dataout2);
