function test_bug2502

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_freqstatistics
% TEST ft_selectdata

% see bug 2502 for more info

% make some data
data1.someparameter = randn(100,10,15);
data1.dimord = 'subj_chan_freq';
data1.freq   = 1:15;
for k = 1:10
  data1.label{k} = ['chan',num2str(k)];
end
data2 = data1;
data2.someparameter = randn(100,10,15);

cfg = [];
cfg.method = 'montecarlo';
cfg.parameter = 'someparameter';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.numrandomization = 0;
cfg.design = [ones(1,100) ones(1,100)*2];
stat = ft_freqstatistics(cfg, data1, data2);

