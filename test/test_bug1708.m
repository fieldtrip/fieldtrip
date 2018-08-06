function test_bug1708

% MEM 1500mb
% WALLTIME 00:10:00

% TEST: test_bug1708
% TEST ft_denoise_synthetic

% reported bug is that ft_denoise_synthetic leads to nans in coilpos and
% coilori

% try to reproduce first
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275.mat'));

cfg  = [];
cfg.gradient = 'none';
data = ft_denoise_synthetic(cfg, data);

cfg = [];
cfg.gradient = 'G3BR';
data = ft_denoise_synthetic(cfg, data);

assert(all(isfinite(data.grad.coilori(:))));
assert(all(isfinite(data.grad.coilpos(:))));
assert(all(isfinite(data.grad.chanori(:))));
assert(all(isfinite(data.grad.chanpos(:))));

%%%% SOMEHOW the test data got lost along the way, therefore it seems
%%%% wisest to just uncomment the next section.
% load('test_bug1708.mat');
% 
% avgData2 = ft_denoise_synthetic(cfg, avgData); % this confirmed the bug
% 
% % however avgData had G1BR balancing, could it be caused by the fact that
% % this balancing is not undone?
% 
% cfg = [];
% cfg.gradient = 'none';
% avgData3 = ft_denoise_synthetic(cfg, avgData);
% cfg.gradient = 'G3BR';
% avgData4 = ft_denoise_synthetic(cfg, avgData3);
% 
% % this indeed seems to be the case
