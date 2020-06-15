function test_bug3461

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY ft_sourcestatistics ft_sourcegrandaverage

%%

% this contains allSourceDiff050
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3461.mat'));

%% Grand average

cfg = [];
cfg.parameter = 'coh';
cfg.keepindividual = 'yes';
avgSourceDiff050 = ft_sourcegrandaverage(cfg, allSourceDiff050{:});

%% Make null condition

cfg = [];
cfg.parameter = 'coh';
cfg.operation = 'multiply';
cfg.scalar = 0;
avgSourceZero = ft_math(cfg, avgSourceDiff050);

%% Run statistics

cfg = [];
cfg.parameter   = 'coh';
cfg.method      = 'montecarlo';
cfg.statistic   = 'depsamplesT';
cfg.alpha       = 0.05;
cfg.tail        = 0;
cfg.correcttail = 'alpha';
cfg.correctm    = 'cluster';
cfg.numrandomization = 1000;

Nsub = 37;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

tval_diff_050 = ft_sourcestatistics(cfg, avgSourceDiff050, avgSourceZero);
