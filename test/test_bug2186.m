function test_bug2186

% MEM 2gb
% WALLTIME 00:10:00

data = [];
for i=1:10
  data.label{i} = num2str(i);
end
for i=1:20
  data.time{i} = (1:1000)/1000;
  data.trial{i} = randn(10, 1000);
end

cfg = [];
cfg.trials = 1:10;
data1 = ft_selectdata(cfg, data);
assert(data1.trial{1}(1)==data.trial{1}(1));


cfg = [];
cfg.trials = 11:20;
data2 = ft_selectdata(cfg, data);
assert(data2.trial{1}(1)==data.trial{11}(1));

%%

cfg = [];
cfg.keeptrials = 'no';
cfg.covariance = 'yes';
timelock1 = ft_timelockanalysis(cfg, data1);
timelock2 = ft_timelockanalysis(cfg, data2);

append12 = ft_appendtimelock([], timelock1, timelock2);
assert(isfield(append12, 'trial')); % renamed from avg
assert(isfield(append12, 'dof'));
assert(isfield(append12, 'var'));
assert(isfield(append12, 'cov')); % this one is new

%%

timelock1a = timelock1;
for i=1:10
  timelock1a.label{i} = num2str(i+10);
end
append1a = ft_appendtimelock([], timelock1, timelock1a);
assert(isequal(size(append1a.cov), [20 20]));
assert(isnan(append1a.cov(20, 1))); % it should be a block diagonal matrix with nan off-diagonal


%%

cfg = [];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes';
timelock1 = ft_timelockanalysis(cfg, data1);
timelock2 = ft_timelockanalysis(cfg, data2);

append12 = ft_appendtimelock([], timelock1, timelock2);
assert(isfield(append12, 'trial')); % avg, dof and var are not present
assert(isfield(append12, 'cov')); % this one is new

%%

timelock1a = timelock1;
for i=1:10
  timelock1a.label{i} = num2str(i+10);
end
append1a = ft_appendtimelock([], timelock1, timelock1a);
assert(isequal(size(append1a.cov), [10 20 20]));
assert(isnan(append1a.cov(1, 20, 1))); % it should be a block diagonal matrix with nan off-diagonal

