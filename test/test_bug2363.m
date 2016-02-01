function test_bug2363

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug2363 ft_selectdata ft_selectdata_new ft_freqstatistics ft_timelockstatistics

%% make some data

data = [];
data.label = {'a','b','c'}'; % note the transpose, this causes the bug
data.dimord = 'chan_freq';
data.freq = 10;
data.powspctrm = randn(numel(data.label),1);

%% select it

cfg = [];
cfg.channel = {'c'};
dat2 = ft_selectdata(cfg, data);

end
