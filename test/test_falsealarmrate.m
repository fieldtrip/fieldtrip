function test_falsealarmrate

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_statistics_montecarlo ft_statistics_analytic ft_statfun_ft_statfun_depsamplesT

% This test script relates to this https://github.com/fieldtrip/website/issues/195

nchan = 100;
ntime = 100;
nsubj = 10; % it failed for 10, it worked for 20

grandavg = [];
grandavg.label = cellstr(num2str((1:nchan)'));
grandavg.time = 1:ntime;
grandavg.avg = randn(2*nsubj, nchan, ntime);
grandavg.dimord = 'subj_chan_time';

%%

cfg = [];
cfg.parameter = 'avg';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.alpha = 0.05;
cfg.correctm = 'no';
cfg.correcttail = 'prob';
cfg.tail = 0;
cfg.numrandomization = 1000;

cfg.design(1,1:2*nsubj) = [ones(1,nsubj) 2*ones(1,nsubj)];
cfg.ivar = 1; % the 1st row in cfg.design contains the independent variable

cfg.method = 'montecarlo';
montecarlo = ft_timelockstatistics(cfg, grandavg);

cfg.method = 'analytic';
analytic = ft_timelockstatistics(cfg, grandavg);

% the false alarm rate should be 5%

assert(mean(montecarlo.mask(:))>0.04);
assert(mean(montecarlo.mask(:))<0.06);

assert(mean(analytic.mask(:))>0.04);
assert(mean(analytic.mask(:))<0.06);

%%

cfg.statistic = 'ft_statfun_depsamplesT';
cfg.design(2,1:2*nsubj) = [1:nsubj 1:nsubj];
cfg.uvar = 2; % the 2nd row in cfg.design contains the subject number

cfg.method = 'montecarlo';
montecarlo = ft_timelockstatistics(cfg, grandavg);

cfg.method = 'analytic';
analytic = ft_timelockstatistics(cfg, grandavg);

% the false alarm rate should be 5%

assert(mean(montecarlo.mask(:))>0.04);
assert(mean(montecarlo.mask(:))<0.06);

assert(mean(analytic.mask(:))>0.04);
assert(mean(analytic.mask(:))<0.06);
