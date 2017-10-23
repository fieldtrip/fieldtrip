function test_bug3229

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_apply_transform ft_componentanalysis ft_rejectcomponent

elec = [];
elec.label   = {'1';'2';'3';'4'};
elec.elecpos = [1 1 1; 2 2 2; 3 3 3; 4 4 4];
elec.chanpos = elec.elecpos;
elec.tra     = eye(4);

%%

bipolar.labelold  = {'1', '2', '3', '4'};
bipolar.labelnew  = {'1', '2', '3'};
bipolar.tra       = [
  +1 -1  0  0
  0 +1 -1  0
  0  0 +1 -1
  ];

elec_bi = ft_apply_montage(elec, bipolar);

% channel names are the same, so keep the same position
assert(isequal(elec_bi.chanpos(1:3,:), elec.chanpos(1:3,:)));

%%

bipolar.labelold  = {'1',   '2',   '3',   '4'};
bipolar.labelnew  = {'1-2', '2-3', '3-4'};
bipolar.tra       = [
  +1 -1  0  0
  0 +1 -1  0
  0  0 +1 -1
  ];

elec_bi = ft_apply_montage(elec, bipolar);

% channel names are not same, so do not keep the same position
assert(~isequal(elec_bi.chanpos(1:3,:), elec.chanpos(1:3,:)));
% channel positions should be half-way
assert(isequal(elec_bi.chanpos(1:3,:), (elec.chanpos(1:3,:)+elec.chanpos(2:4,:))/2));

%%

data = [];
data.label = {'1', '2', '3', '4'};
for i=1:5
  data.trial{i} = randn(4,1000);
  data.time{i}  = (1:1000)/1000;
end
data.elec = elec;
data = ft_checkdata(data);

%%

cfg = [];
cfg.method = 'pca';
cfg.updatesens = 'no';
comp = ft_componentanalysis(cfg, data);
assert(isequal(comp.elec, data.elec)); % should be the same

cfg = [];
cfg.component = []; % keep all
cfg.updatesens = 'no';
backproject = ft_rejectcomponent(cfg, comp);
assert(isequal(backproject.elec, data.elec)); % should still be the same


%%

cfg = [];
cfg.method = 'pca';
cfg.updatesens = 'yes';
comp = ft_componentanalysis(cfg, data);
assert(~isequal(comp.elec, data.elec));

cfg = [];
cfg.component = []; % keep all
cfg.updatesens = 'yes';
backproject = ft_rejectcomponent(cfg, comp);

assert(strcmp(backproject.elec.balance.current, 'invcomp'));
assert(isalmostequal(backproject.elec.tra, eye(4), 'abstol', 1e-9));

cfg = [];
cfg.component = []; % keep all
cfg.updatesens = 'yes';
cleaned = ft_rejectcomponent(cfg, comp, data);

assert(strcmp(cleaned.elec.balance.current, 'reject'));
assert(isalmostequal(cleaned.elec.tra, eye(4), 'abstol', 1e-9));

