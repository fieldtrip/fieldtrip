function test_bug1708

% MEM 1500mb
% WALLTIME 00:03:02

% TEST: test_bug1708
% TEST ft_denoise_synthetic

% reported bug is that ft_denoise_synthetic leads to nans in coilpos and
% coilori

% try to reproduce first
load('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275.mat');

cfg  = [];
cfg.gradient = 'G3BR';
data = ft_denoise_synthetic(cfg, data);

assert(all(isfinite(data.grad.coilori(:))));
assert(all(isfinite(data.grad.coilpos(:))));
assert(all(isfinite(data.grad.chanori(:))));
assert(all(isfinite(data.grad.chanpos(:))));
