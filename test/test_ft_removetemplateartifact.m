function test_ft_removetemplateartifact


% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_removetemplateartifact

fs = 500;
nchan = 32;
start_time = -1; %s
end_time = 2.5; %s
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.time{1} = linspace(start_time, end_time, nsamples);
data.trial{1} = randn(nchan,nsamples);
data.label = cellstr(num2str((1:nchan).'));
data.sens.label = data.label;
data.sens.chanpos = randn(nchan,3);
data.sens.elecpos = data.sens.chanpos;
data.sens.tra = eye(nchan);

template = [];
nsamples = 100;
template.time{1} = linspace(1, nsamples/fs, nsamples);
template.label = data.label;
template.avg = data.trial{1}(:,fs:fs+nsamples-1);
template.var = randn(nchan,nsamples);
template.dimord = 'chan_time';

cfg = [];
cfg.artifact = [fs fs+nsamples-1];
dataout = ft_removetemplateartifact(cfg,data,template);