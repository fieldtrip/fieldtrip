function test_bug1481

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_componentanalysis ft_rejectcomponent ft_apply_montage

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/eeg/preproc_brainvision.mat'));

elec = ft_read_sens('standard_1020.elc');
data.elec = elec;

cfg = [];
cfg.reref = 'yes';
cfg.refchannel = 'all';
data_reref = ft_preprocessing(cfg,data);

cfg = [];
cfg.method = 'fastica';
cfg.numcomponent = 10; % to make it go fast
cfg.randomseed = 13; % so we get the same output each time
comp = ft_componentanalysis(cfg,data_reref);

% this step does not add a balancing matrix to the elec-description
cfg = [];
cfg.component = 2; % chosen randomly
rej1 = ft_rejectcomponent(cfg, comp, data_reref);

% create a montage rereference and call cfg.montage, then assess if
% ft_componentanalysis/ft_rejectcomponent still works correctly.
montage = [];
montage.tra = eye(numel(data.label))-ones(numel(data.label))./numel(data.label);
montage.labelold = data.label;
montage.labelnew = data.label;

cfg = [];
cfg.montage = montage;
data_reref2 = ft_preprocessing(cfg,data);


