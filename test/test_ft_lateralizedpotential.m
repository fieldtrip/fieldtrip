function test_ft_lateralizedpotential

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_lateralizedpotential
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
data.sens.label = data.label;
data.sens.chanpos = randn(nchan,3);
data.sens.elecpos = data.sens.chanpos;
data.sens.tra = eye(nchan);
data.sens.unit = 'cm';

cfg = [];
cfg.channel = data.label(1:nchan/2);
timelockL = ft_timelockanalysis(cfg,data);
cfg.channel = data.label(nchan/2+1:end);
timelockR = ft_timelockanalysis(cfg,data);
ft_checkdata(timelockL, 'datatype', 'timelock');
cfg = [];
cfg.channelcmb = [data.label(1:nchan/2),data.label(nchan/2+1:end)];
dataout = ft_lateralizedpotential(cfg, timelockL, timelockR);

