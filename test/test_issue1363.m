function test_issue1363

% MEM 6gb
% WALLTIME 01:00:00
% DEPENDENCY ft_componentanalysis ft_rejectcomponent

% this test function tests the functionality of the new option cfg.split in
% ft_componentanalysis and ft_rejectcomponent

datadir  = '/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/';
filename = 'sub-15_ses-meg_task-facerecognition_run-01_meg.fif';
dataset  = dccnpath(fullfile(datadir, filename));

cfg         = [];
cfg.dataset = dataset;
cfg.trl     = [60001 110000 0];
cfg.channel = {'MEG';'EEG'};
data        = ft_preprocessing(cfg);

cfg         = [];
cfg.method  = 'pca';
cfg.split   = 'no';
comp1       = ft_componentanalysis(cfg, data);

cfg.split   = 'all';
comp2       = ft_componentanalysis(cfg, data);

cfg.split   = {'MEG' 'EEG'};
comp3       = ft_componentanalysis(cfg, data);

