function test_ft_denoise_dssp

% TEST blabla



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
