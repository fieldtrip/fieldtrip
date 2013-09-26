function test_bug1887

% WALLTIME 00:03:02

% TEST test_bug1887
% TEST ft_checkdata ft_datatype_raw ft_datatype_comp ft_datatype_timelock ft_componentanalysis ft_connectivityanalysis

% this contains raw data, 10 trials with nans
load(dccnfilename('/home/common/matlab/fieldtrip/data/test/bug1887.mat'));

for i=1:10
data.trial{i} = randn(size(data.trial{i}));
end

cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);

cfg = [];
timelock = ft_timelockanalysis(cfg, data);

% simple checks
tmp = ft_checkdata(data, 'datatype', 'raw');
tmp = ft_checkdata(comp, 'datatype', 'comp');
tmp = ft_checkdata(timelock, 'datatype', 'timelock');

% this should also work
tmp = ft_checkdata(comp, 'datatype', 'raw');
assert(isfield(tmp, 'topo'));

% this should also work, but the output should not be comp
tmp = ft_checkdata(comp, 'datatype', 'timelock');
assert(~isfield(tmp, 'topo'));





