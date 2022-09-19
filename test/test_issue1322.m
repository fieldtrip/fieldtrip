function test_issue1322

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_timelocksimulation ft_freqsimulation

headmodel = [];
headmodel.unit = 'm';
headmodel.r = 0.10;
headmodel.o = [0 0 0];
headmodel.cond = 1;

elec = [];
elec.elecpos = randn(10,3);
for i=1:10
  elec.label{i} = num2str(i);
end

%%

cfg = [];
timelock1 = ft_timelocksimulation(cfg);
freq1     = ft_freqsimulation(cfg);

cfg.headmodel = headmodel;
cfg.elec = elec;
dipole1 = ft_dipolesimulation(cfg);

% these do not have to be the same due to different defaults
% assert(isequal(timelock1.time, freq1.time));
% assert(isequal(timelock1.time, dipole1.time));

%%
% specified using cfg.numtrl and cfg.trllen

cfg = [];
cfg.fsample    = 123;
cfg.trllen     = 2;
cfg.numtrl     = 3;
cfg.baseline   = 0.5;
timelock2 = ft_timelocksimulation(cfg);
freq2     = ft_freqsimulation(cfg);

cfg.headmodel = headmodel;
cfg.elec = elec;
dipole2 = ft_dipolesimulation(cfg);

assert(isequal(timelock2.time, freq2.time));
assert(isequal(timelock2.time, dipole2.time));

%%
% specified using cfg.time

cfg = [];
cfg.time = {(0:999)/1000, (0:999)/1000};
timelock3 = ft_timelocksimulation(cfg);
freq3 = ft_freqsimulation(cfg);

cfg.headmodel = headmodel;
cfg.elec = elec;
dipole3 = ft_dipolesimulation(cfg);

assert(isequal(timelock3.time, freq3.time));
assert(isequal(timelock3.time, dipole3.time));

%%
% variable trial lenghths

cfg = [];
cfg.time = {(-100:400)/1000, (-100:900)/1000, (-100:1400)/1000};
timelock4 = ft_timelocksimulation(cfg);
freq4     = ft_freqsimulation(cfg);

cfg.headmodel = headmodel;
cfg.elec = elec;
dipole4 = ft_dipolesimulation(cfg);

assert(isequal(timelock4.time, freq4.time));
assert(isequal(timelock4.time, dipole4.time));

%%

cfg = [];
ft_databrowser(cfg, timelock4);
ft_databrowser(cfg, freq4);
ft_databrowser(cfg, dipole4);

