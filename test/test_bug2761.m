
function test_bug2761

% WALLTIME 00:10:00
% MEM 1GB

% TEST ft_connectivityanalysis ft_connectivity_corr

data = [];
for i=1:5
  data.label{i} = num2str(i);
end

for i=1:13
  data.trial{i} = randn(5,300);
  data.time{i}  = (1:300)/300;
end

cfg = [];
cfg.covariance = 'yes';
timelock = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.method = 'corr';
connectivity = ft_connectivityanalysis(cfg, timelock);

% this failed at the time of reporting this bug
assert(~isfield(connectivity, 'time'));
