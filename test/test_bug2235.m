function test_bug2235

% MEM 500mb
% WALLTIME 00:10:00

% TEST test_bug2235
% TEST ft_denoise_synthetic

fname = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2235');
load(fname);

cfg          = [];
cfg.gradient = 'G3BR';
da_data      = ft_denoise_synthetic(cfg, a_data);