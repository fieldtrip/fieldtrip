function test_pull769

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_read_event

% the data was generated following the example at https://mne.tools/0.15/auto_tutorials/plot_creating_data_structures.html#tut-creating-data-structures
% and written to disk following https://mne.tools/0.15/auto_tutorials/plot_object_epochs.html
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/pull769-epo.fif');

%%

hdr = ft_read_header(filename);
dat = ft_read_data(filename);
evt = ft_read_event(filename);

%%

cfg = [];
cfg.dataset = filename;
cfg.trialdef.eventvalue = 'smiling';
cfg = ft_definetrial(cfg);
data_smiling = ft_preprocessing(cfg);

assert(numel(data_smiling.time)==5);
assert(length(data_smiling.time{1})==200);

cfg = [];
cfg.dataset = filename;
cfg.trialdef.eventvalue = 'frowning';
cfg = ft_definetrial(cfg);
data_frowning = ft_preprocessing(cfg);

assert(numel(data_frowning.time)==5);
assert(length(data_frowning.time{1})==200);
