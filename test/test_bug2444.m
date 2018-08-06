function test_bug2444

% MEM 1000mb
% WALLTIME 00:10:00

% TEST ft_selectdata ft_selectdata_new

%% make some data

comp = [];
comp.label = {'comp1', 'comp2'};
comp.topolabel = {'chan1', 'chan2'};
comp.time = {1:1000, 1:1000};
comp.trial = {randn(2,1000), randn(2,1000)};
comp.topo = eye(2);
comp.unmixing = eye(2);

%% select some data

cfg = [];
cfg.latency = [10 20];
selected = ft_selectdata(cfg, comp);

end
