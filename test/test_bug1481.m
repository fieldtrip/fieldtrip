function test_bug1481

% WALLTIME 00:03:05

% TEST test_bug1481
% TEST ft_rejectcomponent

load('/home/common/matlab/fieldtrip/data/test/latest/raw/eeg/preproc_brainvision.mat');

elec=ft_read_sens('standard_1020.elc');
data.elec=elec;

cfg=[];
cfg.reref='yes';
cfg.refchannel='all';
data_reref=ft_preprocessing(cfg,data);

% TO DO: create a montage rereference and call cfg.montage, then assess if
% ft_rejectcomponent still works correctly.

% cfg=[];
% cfg.montage=montage;
% data_reref=ft_preprocessing(cfg,data);

cfg=[];
cfg.method='fastica';
cfg.numcomponent = 10; % to make it go fast
cfg.randomseed = 13; % so we get the same output each time
comp = ft_componentanalysis(cfg,data_reref);


cfg = [];
cfg.component = 2; % chosen randomly
rej1 = ft_rejectcomponent(cfg, comp, data_reref);


