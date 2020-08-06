function test_bug2855

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_statistics_montecarlo clusterstat

% this is a bug reported long ago (the example data that was created back
% then has never made it here), so it's not clear whether still valid.

% It's about something going wrong when using 'wcm' or 'maxsize', when
% clustering, so it shouldn't depend on the input data type

freq = [];
freq.label = {'a';'b';'c';'d'};
freq.time  = 1:5;
freq.freq  = 1:10;
freq.powspctrm(1:100,1:4,1:10,1:5) = rand(100,4,10,5);
freq.powspctrm(1:50,1:2,1:5,1:3)   = freq.powspctrm(1:50,1:2,1:5,1:3)+0.2;
freq.powspctrm(51:100,3:4,6:10,3:5) = freq.powspctrm(51:100,3:4,6:10,3:5)+0.2;
freq.dimord = 'rpt_chan_freq_time';

n = 50;
cfg=[];
cfg.method = 'montecarlo';
cfg.design = [ones(1,n) ones(1,n)*2];
cfg.ivar = 1;
cfg.numrandomization = 1000;
cfg.correctm = 'no';
cfg.statistic = 'indepsamplesT';
cfg.parameter = 'powspctrm';
stat = ft_freqstatistics(cfg, freq);

cfg.correctm = 'cluster';
for k = 1:4
  cfg.neighbours(k).label = freq.label{k};
  cfg.neighbours(k).neighblabel = {};
end
stat = ft_freqstatistics(cfg, freq);

cfg.clusterstatistic = 'maxsize';
stat = ft_freqstatistics(cfg, freq);

cfg.clusterstatistic = 'wcm';
stat = ft_freqstatistics(cfg, freq);

cfg.clusterstatistic = 'wcm';
cfg.clusterthreshold = 'nonparametric_individual';
stat = ft_freqstatistics(cfg, freq);

cfg.clusterstatistic = 'maxsize';
cfg.clusterthreshold = 'nonparametric_individual';
stat = ft_freqstatistics(cfg, freq);

