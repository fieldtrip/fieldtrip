function test_ft_denoise_pca


% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_denoise_pca

fs = 500;
nchan = 32;
start_time = -1; %s
end_time = 2.5; %s
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.time{1} = linspace(start_time, end_time, nsamples);
data.trial{1} = randn(nchan,nsamples);
data.label = cellstr(num2str((1:nchan).'));

cfg = [];
cfg.refchannel = 'all';
cfg.channel = 'all';
dataout = ft_denoise_pca(cfg,data);
