function test_ft_connectivity_corr

% WALLTIME 00:20:00
% MEM 3gb

tmp.powspctrm = rand(10,3,1);
tmp.freq = 1;
tmp.label = {'1';'2';'3'};
tmp.dimord = 'rpt_chan_freq';
tmp.cumtapcnt = ones(10,1);

cfg = [];
cfg.method = 'powcorr';
stat = ft_connectivityanalysis(cfg, tmp);

isequal(stat.powcorrspctrm,corr(tmp.powspctrm))  
clear tmp stat
    
tmp.trial = rand(10,3,5);
tmp.label = {'1';'2';'3'};
tmp.dimord = 'rpt_chan_time';
tmp.time = [1:5];
tmp.cov = rand(3,3);

cfg = [];
cfg.method = 'corr';
stat = ft_connectivityanalysis(cfg, tmp);
clear tmp stat
