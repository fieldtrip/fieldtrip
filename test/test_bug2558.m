function test_bug2558

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_timelockstatistics

nsubj = 23;
nchan = 274;
ntime = 1560;

gavg_erf7.trial = randn(nsubj,nchan,ntime);
%gavg_erf7.individual = randn(nsubj,nchan,ntime);
gavg_erf7.time = (1:ntime)/1000;
for i=1:nchan
  gavg_erf7.label{i} = num2str(i);
end
gavg_erf7.dimord = 'subj_chan_time';
%gavg_erf7.avg = randn(nchan, ntime);
gavg_erf7.cfg = [];

gavg_erf8.trial = randn(nsubj,nchan,ntime);
%gavg_erf8.individual = randn(nsubj,nchan,ntime);
gavg_erf8.time = (1:ntime)/1000;
for i=1:nchan
  gavg_erf8.label{i} = num2str(i);
end
gavg_erf8.dimord = 'subj_chan_time';
%gavg_erf8.avg = randn(nchan, ntime);
gavg_erf8.cfg = [];

cfg = [];
cfg.channel = 'all'; % modification w.r.t. original bug report
cfg.latency = [0.4 0.7];
cfg.avgovertime = 'no';
cfg.avgoverchan = 'no';
%cfg.avgoverfreq = 'yes'; % <--- is this necessary?
cfg.parameter = 'trial';
cfg.method = 'montecarlo';
%cfg.correctm = 'cluster';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.alpha = 0.05;
%cfg.clusteralpha = 0.05;
%cfg.neighbours  = []; % modification w.r.t. original bug report
cfg.correctm    = 'no';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;

cfg.design(1,1:2*nsubj) =  [ones(1,nsubj) 2*ones(1,nsubj)];
cfg.design(2,1:2*nsubj) = [1:nsubj 1:nsubj];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject numbe

stat = ft_timelockstatistics(cfg, gavg_erf7,gavg_erf8);

