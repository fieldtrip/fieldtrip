function test_rereference

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_preprocessing ft_prepare_montage ft_apply_montage preproc

%% avg

data = [];
data.label = {'1', '2', '3'};
data.time{1} = (1:1000)/1000;
data.trial{1} = rand(length(data.label), length(data.time{1}));

cfg = [];
cfg.implicitref = [];
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = 'all';
data1 = ft_preprocessing(cfg, data);

assert(~isequal(data.trial{1}, data1.trial{1}))
assert(length(data1.label)==3);
assert(all(abs(sum(data1.trial{1},1)) < 10*eps)); % the mean should be zero

cfg = [];
cfg.implicitref = '4';
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = 'all';
data2 = ft_preprocessing(cfg, data);

assert(length(data2.label)==4);
assert(all(abs(sum(data2.trial{1},1)) < 10*eps)); % the mean should be zero

cfg = [];
cfg.implicitref = [];
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = {'1', '2'};
data3 = ft_preprocessing(cfg, data);
assert(all(abs(sum(data3.trial{1}(1:2,:),1)) < 10*eps)); % the mean over the first 2 channels should be zero
assert(~all(abs(mean(data3.trial{1},1)) < 10*eps)); % the mean over all channels should NOT be zero

%% median

data = [];
data.label = {'1', '2', '3'};
data.time{1} = (1:1000)/1000;
data.trial{1} = rand(length(data.label), length(data.time{1}));

cfg = [];
cfg.implicitref = [];
cfg.reref = 'yes';
cfg.refmethod = 'median';
cfg.refchannel = 'all';
data1 = ft_preprocessing(cfg, data);

assert(~isequal(data.trial{1}, data1.trial{1}))
assert(length(data1.label)==3);
assert(all(abs(median(data1.trial{1},1)) < 10*eps)); % the median should be zero
assert(~all(abs(mean(data1.trial{1},1)) < 10*eps)); % the mean should NOT be zero

cfg = [];
cfg.implicitref = '4';
cfg.reref = 'yes';
cfg.refmethod = 'median';
cfg.refchannel = 'all';
data2 = ft_preprocessing(cfg, data);

assert(length(data2.label)==4);
assert(all(abs(median(data2.trial{1},1)) < 10*eps)); % the median should be zero
assert(~all(abs(mean(data2.trial{1},1)) < 10*eps)); % the mean should NOT be zero

cfg = [];
cfg.implicitref = [];
cfg.reref = 'yes';
cfg.refmethod = 'median';
cfg.refchannel = {'1', '2'};
data3 = ft_preprocessing(cfg, data);
assert(all(abs(median(data3.trial{1}(1:2,:),1)) < 10*eps)); % the median over the first 2 channels should be zero
assert(~all(abs(median(data3.trial{1},1)) < 10*eps)); % the median over all channels should NOT be zero
assert(~all(abs(mean(data3.trial{1},1)) < 10*eps)); % the mean over all channels should NOT be zero

%% rest

elec = [];
elec.label = arrayfun(@num2str, 1:162, 'UniformOutput', false)';
elec.elecpos = mesh_sphere(162);

headmodel = [];
headmodel.type = 'singlesphere';
headmodel.r = 1;
headmodel.o = [0 0 0];
headmodel.cond = 1;

cfg = [];
cfg.headmodel = headmodel;
cfg.elec = elec;
cfg.method = 'basedonvol';
cfg.inwardshift = 0.1;
sourcemodel = ft_prepare_sourcemodel(cfg);

if false
  figure
  ft_plot_headmodel(headmodel);
  alpha 0.5
  ft_plot_mesh(sourcemodel);
  ft_plot_sens(elec, 'label', 'label', 'elecshape', 'disc');
end

cfg = [];
cfg.headmodel = headmodel;
cfg.elec = elec;
cfg.sourcemodel = sourcemodel;
leadfield = ft_prepare_leadfield(cfg);

lfmat = cat(2, leadfield.leadfield{:});
assert(isequal(size(lfmat), [162, 642*3]));
assert(~any(isnan(lfmat(:))));
assert(all(abs(mean(lfmat,1)) < 10*eps)); % the leadfield should be average referenced

data = [];
data.label = arrayfun(@num2str, 1:162, 'UniformOutput', false)';
data.time{1} = (1:1000)/1000;
data.trial{1} = rand(length(data.label), length(data.time{1}));

cfg = [];
cfg.implicitref = [];
cfg.reref = 'yes';
cfg.refmethod = 'rest';
cfg.refchannel = 'all';
cfg.leadfield = leadfield;
data1 = ft_preprocessing(cfg, data);

% the REST rereferenced data should have a mean that is close to zero
% since the electrode and dipole distribution are very homogenous

assert(all(abs(mean(data1.trial{1},1)) < 10*eps));

%% bipolar

nelec = 10; % per shaft
nshaft = 5;

shaft = {'A', 'B', 'C', 'D', 'E'};

data = [];
data.label = {};
for i=1:nshaft
  for j=1:nelec
    data.label{end+1} = sprintf('%s%d', shaft{i}, j);
  end
end
data.time{1} = (1:1000)/1000;
data.trial{1} = rand(length(data.label), length(data.time{1}));

% note that the ordering should be 1, 2, 3, ... 10, and not 1, 10, 2, 3, 4, ...

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'bipolar';
cfg.groupchans = 'no';
data1 = ft_preprocessing(cfg, data);
assert(length(data1.label)==49); % 50 channels, hence 49 differences

a1_min_a2 = data.trial{1}(1,:) - data.trial{1}(2,:);
assert(all(abs(data1.trial{1}(1,:)-a1_min_a2)<10*eps));

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'bipolar';
cfg.groupchans = 'yes';
data2 = ft_preprocessing(cfg, data);
assert(length(data2.label)==45); % 5x10 channels, hence 5x9 differences

a1_min_a2 = data.trial{1}(1,:) - data.trial{1}(2,:);
assert(all(abs(data2.trial{1}(1,:)-a1_min_a2)<10*eps));


%% laplace

nelec = 10; % per shaft
nshaft = 5;

shaft = {'A', 'B', 'C', 'D', 'E'};

data = [];
data.label = {};
for i=1:nshaft
  for j=1:nelec
    data.label{end+1} = sprintf('%s%d', shaft{i}, j);
  end
end
data.time{1} = (1:1000)/1000;
data.trial{1} = rand(length(data.label), length(data.time{1}));

% note that the ordering should be 1, 2, 3, ... 10, and not 1, 10, 2, 3, 4, ...

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'laplace';
cfg.groupchans = 'no';
data1 = ft_preprocessing(cfg, data);

assert(length(data1.label)==50);

cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'laplace';
cfg.groupchans = 'yes';
data2 = ft_preprocessing(cfg, data);

assert(length(data1.label)==50);

% channel 1 is in both cases on the edge, hence the same
% channel 10 and 11 are in one case in the middle, in the other at two separate edges of the shafts
% channel 2-9 are always embedded in the shaft, so the same

assert( isequal(data1.trial{1}(1,:), data2.trial{1}(1,:)));
assert(~isequal(data1.trial{1}(10,:), data2.trial{1}(10,:)));
assert(~isequal(data1.trial{1}(11,:), data2.trial{1}(11,:)));
assert( isequal(data1.trial{1}(2:9,:), data2.trial{1}(2:9,:)));
