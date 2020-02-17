function test_issue1294

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_statfun_bayesfactor

%%

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/effectsize/ERF_orig.mat'));

cfg = [];
cfg.keepindividual = 'yes';
grandavgFIC = ft_timelockgrandaverage(cfg, allsubjFIC{:});
grandavgFC  = ft_timelockgrandaverage(cfg, allsubjFC{:});

%%

cfg = [];
cfg.channel = {'MLT12', 'MLT13', 'MLT23', 'MLT24', 'MLT32', 'MLT33', 'MLT41'};
cfg.latency = [0.35 0.55];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.method = 'analytic';
cfg.statistic = 'bayesfactor'; % see FT_STATFUN_BAYESFACTOR
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  ];
effect_avg = ft_timelockstatistics(cfg, grandavgFIC, grandavgFC);


%%

cfg = [];
cfg.parameter = 'individual';
cfg.channel = 'MEG';
cfg.method = 'analytic';
cfg.statistic = 'bayesfactor'; % see FT_STATFUN_BAYESFACTOR
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10
  ];
effect_all = ft_timelockstatistics(cfg, grandavgFIC, grandavgFC);

cfg.statistic = 'depsamplesT';
effect_tstat = ft_timelockstatistics(cfg, grandavgFIC, grandavgFC);


%% 
% this one is only for testing purposes
% the data should be paired for a proper analysis

cfg = [];
cfg.parameter = 'individual';
cfg.channel = 'MEG';
cfg.method = 'analytic';
cfg.statistic = 'bayesfactor'; % see FT_STATFUN_BAYESFACTOR
cfg.latency = [0 0.4];
cfg.ivar = 1;
cfg.design = [
  1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
  ];
effect_nonpaired = ft_timelockstatistics(cfg, grandavgFIC, grandavgFC);


%%

cfg = [];
cfg.layout = 'CTF151_helmet.mat';
cfg.parameter = 'bf10';
cfg.ylim = [0 10];
figure; ft_multiplotER(cfg, effect_all); title('Bayes Factor')

%%

cfg = [];
cfg.layout = 'CTF151_helmet.mat';
cfg.parameter = 'tstat';
cfg.ylim = [-6 6];
figure; ft_multiplotER(cfg, effect_all); title('BF tstat')

cfg.parameter = 'stat';
figure; ft_multiplotER(cfg, effect_tstat); title('regular stat') 

