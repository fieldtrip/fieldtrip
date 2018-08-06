function test_bug3048

% MEM 2gb
% WALLTIME 00:20:00

% TEST ft_preamble ft_preamble_randomseed ft_dipolesimulation
% TEST ft_freqsimulation ft_connectivitysimulation ft_statistics_montecarlo
% TEST ft_freqstatistics ft_timelockstatistics

%%
% test ft_freqsimulation
cfg = [];

freq0 = ft_freqsimulation(cfg);
freq1 = ft_freqsimulation(cfg);

cfg.randomseed = freq1.cfg.callinfo.randomseed;
freq2 = ft_freqsimulation(cfg);

assert(~isequal(freq0.trial{1}, freq1.trial{1}));
assert( isequal(freq1.trial{1}, freq2.trial{1}));

%%
% test ft_connectivitysimulation
cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';
cfg.params(:,:,1) = [
  0.8  0.0  0.0;
  0.0  0.9  0.5 ;
  0.4  0.0  0.5];
cfg.params(:,:,2) = [
 -0.5  0.0  0.0;
  0.0 -0.8  0.0;
  0.0  0.0 -0.2];
cfg.noisecov = [
  0.3  0.0  0.0;
  0.0  1.0  0.0;
  0.0  0.0  0.2];

conn0 = ft_connectivitysimulation(cfg);
conn1 = ft_connectivitysimulation(cfg);

cfg.randomseed = conn1.cfg.callinfo.randomseed;
conn2 = ft_connectivitysimulation(cfg);

assert(~isequal(conn0.trial{1}, conn1.trial{1}));
assert( isequal(conn1.trial{1}, conn2.trial{1}));

%%
% test ft_dipolesimulation
cfg = [];
cfg.relnoise = 1;
cfg.dip.pos = [0 0 0.05];
cfg.headmodel.unit = 'm';
cfg.headmodel.r = 0.10;
cfg.headmodel.o = [0 0 0];
cfg.headmodel.cond = 1;
cfg.elec.elecpos = randn(10,3);
for i=1:10
  cfg.elec.label{i} = num2str(i);
end

dip0 = ft_dipolesimulation(cfg);
dip1 = ft_dipolesimulation(cfg);

cfg.randomseed = dip1.cfg.callinfo.randomseed;
dip2 = ft_dipolesimulation(cfg);

assert(~isequal(dip0.trial{1}, dip1.trial{1}));
assert( isequal(dip1.trial{1}, dip2.trial{1}));

%%
% test ft_statistics_montecarlo
cfg = [];
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 100;
design = [ones(1,10) ones(1,10)*2;1:10 1:10];
dat    = randn(5,20);

[stat1, cfg1] = ft_statistics_montecarlo(cfg, dat, design);
[stat2, cfg2] = ft_statistics_montecarlo(cfg, dat, design);
cfg.randomseed = cfg2.callinfo.randomseed;
[stat3, cfg3] = ft_statistics_montecarlo(cfg, dat, design);

assert(~isequal(stat1.prob, stat2.prob));
assert( isequal(stat2.prob, stat3.prob));

%% 
% test the higher level stats functions
cfg.randomseed = [];
cfg.design     = design;
cfg.method     = 'montecarlo';
freq.powspctrm = randn(20,5,10);
freq.freq      = 1:10;
freq.dimord    = 'rpt_chan_freq';
for k = 1:5, freq.label{k} = sprintf('chan%02d',k); end
stat1 = ft_freqstatistics(cfg, freq);
stat2 = ft_freqstatistics(cfg, freq);
cfg.randomseed = stat2.cfg.callinfo.randomseed;
stat3 = ft_freqstatistics(cfg, freq);
assert(~isequal(stat1.prob, stat2.prob));
assert( isequal(stat2.prob, stat3.prob));

cfg.randomseed = [];
tlck.trial     = randn(20,5,10);
tlck.time      = 1:10;
tlck.label     = freq.label;
tlck.dimord    = 'rpt_chan_time';
stat1 = ft_timelockstatistics(cfg, tlck);
stat2 = ft_timelockstatistics(cfg, tlck);
cfg.randomseed = stat2.cfg.callinfo.randomseed;
stat3 = ft_timelockstatistics(cfg, tlck);
assert(~isequal(stat1.prob, stat2.prob));
assert( isequal(stat2.prob, stat3.prob));

