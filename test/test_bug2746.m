function test_bug2746

% TEST test_bug2746
% TEST loadcnt

% MEM 150mb
% WALLTIME 00:10:00

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/neuroscan32/Subject1_MP.cnt');

cfg.dataset = filename;
cfg.trl     = [1 697 0;15335 16031 0;16032 16728 0];
data = ft_preprocessing(cfg);

assert(~isequal(data.trial{1},data.trial{3}));
