function test_bug2585

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_componentanalysis ft_preamble_randomseed

p = 5;
n = 100;
data = [];
for i=1:p
  data.label{i} = num2str(i);
end
data.trial = {randn(p, n).^2}; % square it, it should be non-normal
data.time = {1:n};

% perform the independent component analysis (i.e., decompose the data)
cfg = [];
cfg.channel = {'all'};

% cfg.randomseed = 1;

cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
% cfg.numcomponent = 20;
% cfg.runica.pca = 20;

comp1 = ft_componentanalysis(cfg, data);

cfg.randomseed = comp1.cfg.callinfo.randomseed;
comp2 = ft_componentanalysis(cfg, data);

cfg.randomseed = comp1.cfg.callinfo.randomseed;
comp3 = ft_componentanalysis(cfg, data);

cfg.randomseed = 42;
comp4 = ft_componentanalysis(cfg, data);

assert( isequal(comp1.trial{1}, comp2.trial{1})); % this test fails unexpectedly, which explains the bug report
assert( isequal(comp1.trial{1}, comp3.trial{1}));
assert(~isequal(comp1.trial{1}, comp4.trial{1}));

