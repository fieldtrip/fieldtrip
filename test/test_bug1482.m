% function test_bug1482

% WALLTIME 00:10:00
% MEM1gb

elec = [];
elec.label = {'1', '2', '3', '4'};
elec.unit = 'cm';
elec.elecpos = [
  1 1 1
  2 2 2
  3 3 3
  4 4 4
  ];

data = [];
data.label = {'1', '2', '3', '4', '5'}; % one more than elec
% data.label = {'1', '2', '3', '4'};
for i=1:10
  data.time{i} = (1:1000)/1000;
  data.trial{i} = randn(length(data.label), 1000);
end
data.elec = elec;

%%

tmpcfg = [];
tmpcfg.reref = 'yes';
tmpcfg.implicitref = [];
tmpcfg.refchannel = 'all';

cfg = [];
cfg.montage = ft_prepare_montage(tmpcfg, data);
cfg.updatesens = 'yes';
data_reref1 = ft_preprocessing(cfg, data);

%%

cfg = [];
cfg.reref = 'yes';
cfg.implicitref = [];
cfg.refchannel = 'all';
cfg.updatesens = 'yes';
data_reref2 = ft_preprocessing(cfg, data);

%%

assert(isalmostequal(data_reref1.trial{1}, data_reref2.trial{1}, 'reltol', 1e-6));
assert(isalmostequal(data_reref1.elec, data_reref2.elec, 'reltol', 1e-6));

%%

tmpcfg = [];
tmpcfg.reref = 'yes';
tmpcfg.implicitref = 'REF';
tmpcfg.refchannel = 'all';

cfg = [];
cfg.updatesens = 'yes';
cfg.montage = ft_prepare_montage(tmpcfg, data);
data_reref1 = ft_preprocessing(cfg, data);

%%

cfg = [];
cfg.reref = 'yes';
cfg.implicitref = 'REF';
cfg.refchannel = 'all';
cfg.updatesens = 'yes';
data_reref2 = ft_preprocessing(cfg, data);

%%

assert(isalmostequal(data_reref1.trial{1}, data_reref2.trial{1}, 'reltol', 1e-6));
assert(isalmostequal(data_reref1.elec, data_reref2.elec, 'reltol', 1e-6));



