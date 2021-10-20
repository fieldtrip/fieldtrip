function test_ft_appenddata

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_appenddata

fs = 500;
nchan = 32;
start_time = -1; %s
end_time = 2.5; %s
nsamples = (end_time - start_time) * fs + 1;

data1 = [];
data1.time{1} = linspace(start_time, end_time, nsamples);
data1.trial{1} = randn(nchan,nsamples);
data1.sampleinfo = [1 nsamples];
data1.label = cellstr(num2str((1:nchan).'));
data1.sens.label = data1.label;
data1.sens.chanpos = randn(nchan,3);
data1.sens.elecpos = data1.sens.chanpos;
data1.sens.tra = eye(nchan);

% test case of different channels, same time frames
data2 = data1;
data2.trial{1} = randn(nchan,nsamples);
data2.label = cellstr(num2str((nchan+1:2*nchan).'));
data2.sens.label = data2.label;
data2.sens.chanpos = randn(nchan,nsamples);
data2.sens.elecpos = data2.sens.chanpos;
cfg = [];
dataout = ft_appenddata(cfg,data1,data2);

% test case of different time frames, same channels
data2 = data1;
nsamples = fs-1;
data2.time{1} = linspace(end_time, end_time+1, nsamples);
data2.trial{1} = randn(nchan,nsamples);
data2.sampleinfo = [1 nsamples];
cfg = [];
dataout = ft_appenddata(cfg,data1,data2);
