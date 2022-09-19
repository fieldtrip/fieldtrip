function test_ft_channelnormalise

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_channelnormalise

fs = 500;
nchan = 32;
start_time = -1; % seconds
end_time = 2.5; % seconds
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.time{1} = linspace(start_time, end_time, nsamples);
data.trial{1} = randn(nchan,nsamples);
data.label = cellstr(num2str((1:nchan).'));

cfg = [];
dataout = ft_channelnormalise(cfg, data);
