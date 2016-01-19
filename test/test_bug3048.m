function test_bug3048

% MEM 2gb
% WALLTIME 00:20:00

% TEST ft_preamble ft_preamble_randomseed ft_dipolesimulation ft_freqsimulation ft_connectivitysimulation


%%

cfg = [];

freq0 = ft_freqsimulation(cfg);
freq1 = ft_freqsimulation(cfg);

cfg.randomseed = freq1.cfg.callinfo.randomseed;
freq2 = ft_freqsimulation(cfg);

assert(~isequal(freq0.trial{1}, freq1.trial{1}));
assert( isequal(freq1.trial{1}, freq2.trial{1}));

%%

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



