function test_ft_prepare_neighbours

% TEST test_ft_prepare_neighbours
% TEST ft_prepare_neighbours

datainfo = ref_datasets;

% get an MEG and an EEG set (hard-coded
eeginfo = datainfo(3);
meginfo = datainfo(7);

% do the MEG processing
fname = fullfile(meginfo.origdir,'latest', 'raw',meginfo.type,['preproc_' meginfo.datatype]);
load(fname);

%% triangulation method
cfg = [];
cfg.method = 'triangulation';
neighbours = ft_neighbourselection(cfg, data);

% do only on some channels
cfg.channel = data.label(1:end/2);
neighbours = ft_neighbourselection(cfg, data);

% remove some channels
data.label(1) = [];
neighbours = ft_neighbourselection(cfg, data);

%% template method
cfg = [];
cfg.method = 'template';
tic
neighbours = ft_neighbourselection(cfg, data);
toc

cfg = [];
cfg.method = 'template';
cfg.template = 'bti248_neighb.mat';
tic
neighbours = ft_neighbourselection(cfg);
toc

% do only on some channels
cfg.channel = data.label(1:end/2);
neighbours = ft_neighbourselection(cfg, data);

% remove some channels
data.label(1) = [];
neighbours = ft_neighbourselection(cfg, data);
