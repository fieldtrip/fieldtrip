function test_ft_prepare_montage

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

%%

fsample = 1000;
ntrial = 20;
nsample = 1000;
nchan = 7;

data = [];
for i=1:nchan
  data.label{i} = sprintf('%d', i);
end
for i=1:ntrial
  data.time{i}  = (1:nsample)./fsample;
  data.trial{i} = randn(nchan, nsample);
  data.trial{i}(:,1) = 1:nchan;  % replace the first sample
end

cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);


%%

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = 'all';

data1a = ft_preprocessing(cfg, data);
cfg.implicitref = '0';
data1b = ft_preprocessing(cfg, data);

assert(numel(data1a.label)==nchan);
assert(isalmostequal(mean(data1a.trial{1}), zeros(1,1000), 'abstol', 100*eps));

assert(numel(data1b.label)==nchan+1);
assert(isalmostequal(mean(data1b.trial{1}), zeros(1,1000), 'abstol', 100*eps));

%%

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'median';
cfg.refchannel = 'all';

data2a = ft_preprocessing(cfg, data);
cfg.implicitref = '0';
data2b = ft_preprocessing(cfg, data);

assert(numel(data2a.label)==nchan);
assert(isalmostequal(median(data2a.trial{1}), zeros(1,1000), 'abstol', 100*eps));

assert(numel(data2b.label)==nchan+1);
assert(isalmostequal(median(data2b.trial{1}), zeros(1,1000), 'abstol', 100*eps));

%%

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'bipolar';
data3a = ft_preprocessing(cfg, data);
cfg.implicitref = '0';
data3b = ft_preprocessing(cfg, data);

assert(numel(data3a.label)==nchan-1);
assert(all(data3a.trial{1}(:,1)==-1)) % this is 2-1, 3-2, etc.

assert(numel(data3b.label)==nchan+1-1);
assert(all(data3b.trial{1}(1:nchan-1,1)==-1)) % this is 2-1, 3-2, etc.
assert(all(data3b.trial{1}(nchan,1)==nchan))  % this is nchan-0

%%

cfg = [];
cfg.refmethod = 'comp';
montage = ft_prepare_montage(cfg, comp);

cfg = [];
cfg.montage = montage;
data4a = ft_preprocessing(cfg, data);

cfg = [];
cfg.refmethod = 'invcomp';
montage = ft_prepare_montage(cfg, comp);

cfg = [];
cfg.montage = montage;
data4b = ft_preprocessing(cfg, comp);
% FIXME I would expect the output not to contain topo/unmixing/topolabel

%%

fsample = 1000;
ntrial = 20;
nsample = 1000;

data1020 = [];
data1020.label = ft_senslabel('eeg1020');

nchan = length(data1020.label);

for i=1:ntrial
  data1020.time{i}  = (1:nsample)./fsample;
  data1020.trial{i} = randn(nchan, nsample);
end

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'doublebanana';
data1020_doublebanana = ft_preprocessing(cfg, data1020);

cfg.refmethod = 'longitudinal';
data1020_longitudinal = ft_preprocessing(cfg, data1020);

cfg.refmethod = 'transverse';
data1020_transverse = ft_preprocessing(cfg, data1020);

cfg.refmethod = 'circumferential';
data1020_circumferential = ft_preprocessing(cfg, data1020);

