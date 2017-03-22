function test_bug2553

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_componentanalysis

%% generate some data

data = [];
data.label = {'a' 'b' 'c'}';
data.time = {1:1000 1:1000 1:1000};
data.trial = {randn(3,1000) randn(3,1000) randn(3,1000)};
data.sampleinfo = [1 1000;
  1501 2500;
  3501 4500];
data.trialinfo = (1:3)';

%% componentanalysis

cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);

%% verify

assert(isequal(data.sampleinfo, comp.sampleinfo) && isequal(data.trialinfo, comp.trialinfo));
