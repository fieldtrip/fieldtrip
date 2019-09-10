function test_issue1162

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_componentanalysis ft_rejectcomponent ft_apply_montage

%%

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/raw.fif');
cfg.trl = [1 7500 0]; % 30 seconds @250Hz
cfg.continuous = 'yes';
cfg.channel = {'MEG'};
data = ft_preprocessing(cfg);

cfg = [];
cfg.numcomponent = 12;
cfg.method = 'pca'; % the type of unmixing does not matter, this is the fastest
comp = ft_componentanalysis(cfg, data);

% this originally failed ...
cfg = [];
cfg.component = [1 2 3];
clean = ft_rejectcomponent(cfg, comp, data);

%%

% prevent data from being all zero
data1 = data;
data1.trial{1}(1:3,:) = randn(size(data1.trial{1}(1:3,:)));

cfg = [];
cfg.numcomponent = 12;
cfg.method = 'pca'; % the type of unmixing does not matter, this is the fastest
comp1 = ft_componentanalysis(cfg, data1);

% ... whereas this worked
cfg = [];
cfg.component = [1 2 3];
clean1 = ft_rejectcomponent(cfg, comp1, data1);
