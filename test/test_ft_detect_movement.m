function test_ft_detect_movement


% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_detect_movement

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
cfg.velocity2D.kernel = [1 1 0 -1 -1].*(fs/6);
[cfg,movement] = ft_detect_movement(cfg,data);
