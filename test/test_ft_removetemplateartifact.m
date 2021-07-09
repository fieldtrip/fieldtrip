function test_ft_removetemplateartifact

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_removetemplateartifact

fs = 500;
nchan = 32;
start_time = -1; % seconds
end_time = 2.5; % seconds
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.time{1} = linspace(start_time, end_time, nsamples);
data.trial{1} = randn(nchan,nsamples);
data.label = cellstr(num2str((1:nchan).'));

nsamples = 100;

template = [];
template.time{1} = linspace(1, nsamples/fs, nsamples);
template.label = data.label;
template.avg = data.trial{1}(:,fs:fs+nsamples-1);
template.var = randn(nchan,nsamples);
template.dimord = 'chan_time';

cfg = [];
cfg.artifact = [fs fs+nsamples-1];
dataout = ft_removetemplateartifact(cfg,data,template);
