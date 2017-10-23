function test_bug2593

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_componentanalysis

% this function tests the cfg.numcomponent and related options for
% ft_componentanalysis and its 3rd-party subfunctions
% note: binica is not tested, since it is typically not used (broken for
% some other reason)

%% create some data

data = [];
data.label = {'a','b','c'};
data.trial = {randn(3,1000),randn(3,1000),randn(3,1000)};
data.time = {1:1000 1:1000 1:1000};

%%

cfg = [];
cfg.method = 'runica';
comp = ft_componentanalysis(cfg, data);
assert(numel(comp.label) == 3);

%%

cfg = [];
cfg.method = 'fastica';
comp = ft_componentanalysis(cfg, data);
assert(numel(comp.label) == 3);

%%

cfg = [];
cfg.method = 'runica';
cfg.numcomponent = 2;
comp = ft_componentanalysis(cfg, data);
assert(numel(comp.label) == 2);

%%

cfg = [];
cfg.method = 'fastica';
cfg.numcomponent = 2;
comp = ft_componentanalysis(cfg, data);
assert(numel(comp.label) == 2);

%%

cfg = [];
cfg.method = 'runica';
cfg.runica.pca = 2;
comp = ft_componentanalysis(cfg, data);
assert(numel(comp.label) == 2);

%%

cfg = [];
cfg.method = 'fastica';
cfg.fastica.numOfIC = 2;
comp = ft_componentanalysis(cfg, data);
assert(numel(comp.label) == 2);

%%

cfg = [];
cfg.method = 'runica';
cfg.runica.pca = 2;
cfg.numcomponent = 2;

ok = 0;
try
  comp = ft_componentanalysis(cfg, data);
catch err
  ok = 1;
end

if ~ok
  error('ft_componentanalysis should have thrown an error');
end

%%

cfg = [];
cfg.method = 'fastica';
cfg.fastica.numOfIC = 2;
cfg.numcomponent = 2;

ok = 0;
try
  comp = ft_componentanalysis(cfg, data);
catch err
  ok = 1;
end

if ~ok
  error('ft_componentanalysis should have thrown an error');
end
end
