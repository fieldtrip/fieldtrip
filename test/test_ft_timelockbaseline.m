function test_ft_timelockbaseline

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_timelockbaseline
% DATATYPE timelock

fs = 500;
nchan = 32;
start_time = -1; %s
end_time = 2.5; %s
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.time{1} = linspace(start_time, end_time, nsamples);
data.trial{1} = randn(nchan,nsamples);
data.sampleinfo = [1 nsamples];
data.label = cellstr(num2str((1:nchan).'));

cfg = [];
timelock = ft_timelockanalysis(cfg,data);
ft_checkdata(timelock, 'datatype', 'timelock');

cfg = [];
cfg.baseline = [-1 0];
cfg.channel = 'all';
dataout = ft_timelockbaseline(cfg, timelock);

