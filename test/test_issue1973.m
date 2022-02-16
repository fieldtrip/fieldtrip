function test_issue1973

% WALLTIME 00:10:00
% MEM 4gb
% DEPENDENCY ft_statfun_actvsblT

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1973.mat'));

%%

ntrials = size(tlprestim.trial,1);
design  = zeros(2,2*ntrials);
design(1,1:ntrials) = 2; %NOTE THAT THIS IS SWAPPED WRT CHARLOTTE'S .MLX
design(1,ntrials+1:2*ntrials) = 1;
design(2,1:ntrials) = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];

origtrialdims         = size(tlprestim.trial);
tlprestimhacked       = tlprestim;
tlprestimhacked.trial = repmat(mean(tlprestim.trial,3),[1 1 origtrialdims(3)]);
tlprestimhacked.time  = tlpoststim.time;

% Cluster-based permutation

cfg = [];
%cfg.latency = [0 0.200];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_actvsblT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% prepare_neighbours determines what sensors may form clusters
cfg.neighbours       = neighbours;

cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;


tlprestim.time  = tlpoststim.time; % this is needed, because FieldTrip checks for the overlap in the time axes

% this originally gave an error because of an issue with ft_statfun_actvsblT
[timelockstats1] = ft_timelockstatistics(cfg, tlprestim, tlpoststim);

[timelockstats2] = ft_timelockstatistics(cfg, tlprestimhacked, tlpoststim);

cfg.statistic         = 'ft_statfun_depsamplesT';

[timelockstats3] = ft_timelockstatistics(cfg, tlprestimhacked,tlpoststim);

assert(isalmostequal(timelockstats1.stat,timelockstats2.stat, 'abstol', 10*eps));
assert(isalmostequal(timelockstats1.stat,timelockstats3.stat, 'abstol', 10*eps));


