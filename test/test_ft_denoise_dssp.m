function test_ft_denoise_dssp

% MEM 3000mb
% WALLTIME 00:10:00

% TEST ft_denoise_dssp

% below is based on simulated data using Kensuke's code, and should behave as the
% bdssp_main script that he uses to demonstrate the algorithm
load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull800/data.mat'));
tmp = load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull800/grid.mat'));

cfg = [];
cfg.grid = tmp.grid;
cfg.dssp.n_space = 30;
cfg.dssp.n_in    = 30;
cfg.dssp.n_out   = 30;
dataout = ft_denoise_dssp(cfg, data);

load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull800/data_lowsnr.mat'));
dataout2 = ft_denoise_dssp(cfg, data);
