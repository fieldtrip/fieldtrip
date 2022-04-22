function test_example_threshold_free_cluster_enhancement_20210618

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY ft_statisatics_montecarlo

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/threshold_free_cluster_enhancement/ERF_orig.mat'));

%%

Nsubj       = 10;
design      = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];                % this is the uvar (unit-of-observation variable)
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];  % this is the ivar (independent variable)

% channel neighbours are not really needed for the subsequent example, since we restrict the test to a single channel
% nevertheless, we will determine the channel neighbours here anyway
cfg          = [];
cfg.method   = 'template';
cfg.template = 'CTF151_neighb.mat';
nb           = ft_prepare_neighbours( cfg, allsubjFIC{1} );

%% Conventional cluster-based permutation test
cfg                  = [];
cfg.design           = design;
cfg.uvar             = 1;
cfg.ivar             = 2;
cfg.channel          = {'MLT12'};
cfg.latency          = [0 1];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.numrandomization = 'all'; % there are 10 subjects, so 2^10=1024 possible raondomizations
cfg.tail             = 0;
cfg.alpha            = 0.025; % since we are testing two tails
cfg.neighbours       = nb;
cfg.correctm         = 'cluster';
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_individual';
cfg.clustertail      = 0;
cfg.clusteralpha     = 0.01;
stat01 = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});

cfg.clusteralpha     = 0.05;
stat05 = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});

%% Cluster-based permutation test with TFCE method
cfg                  = [];
cfg.design           = design;
cfg.uvar             = 1;
cfg.ivar             = 2;
cfg.channel          = {'MLT12'};
cfg.latency          = [0 1];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.numrandomization = 'all';  % there are 10 subjects, so 2^10=1024 possible raondomizations
cfg.tail             = 0;
cfg.alpha            = 0.025;
cfg.neighbours       = nb;
cfg.correctm         = 'tfce';
cfg.tfce_H           = 2;      % default setting
cfg.tfce_E           = 0.5;    % default setting
statA = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});

cfg.tfce_H           = 2;
cfg.tfce_E           = 0.25;
statB = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:});

%% Visualize the results
figure(1), clf, hold on,
set(gcf, 'units','centimeters','position',[0 0 18 10] );
subplot(2,2,1), hold on, grid on,
title( 'TFCE: H=2, E=0.5' );
plot( statA.time, statA.stat_tfce, 'r');
plot( statA.time(statA.mask), 400*ones(sum(statA.mask),1), 'r', 'linewidth',2 );
ylabel( 'TFCE' );
ylim( [-20, 420] );
xticks( 0:0.2:1 );
set(gca, 'TickDir','out' );

subplot(2,2,3), hold on, grid on,
title( 'TFCE: H=2, E=0.25' );
plot( statB.time, statB.stat_tfce, 'r');
plot( statB.time(statB.mask), 400*ones(sum(statB.mask),1), 'r', 'linewidth',2 );
ylabel( 'TFCE' );
ylim( [-20, 420] );
xticks( 0:0.2:1 );
set(gca, 'TickDir','out' );

subplot(2,2,2), hold on, grid on,
title( 'clusteralpha = .01' );
plot( stat01.time, stat01.stat, 'k');
plot( stat01.time(stat01.mask), 7.5*ones(sum(stat01.mask),1), 'k', 'linewidth',2 );
ylabel( 't-value' );
ylim( [-3, 8] );
xticks( 0:0.2:1 );
set(gca, 'TickDir','out' );

subplot(2,2,4), hold on, grid on,
title( 'clusteralpha = .05' );
plot( stat05.time, stat05.stat, 'k');
plot( stat05.time(stat05.mask), 7.5*ones(sum(stat05.mask),1), 'k', 'linewidth',2 );
ylabel( 't-value' );
ylim( [-3, 8] );
xticks( 0:0.2:1 );
set(gca, 'TickDir','out' );
