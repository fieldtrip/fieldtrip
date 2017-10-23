function test_bug2746

% TEST loadcnt

% MEM 1500mb
% WALLTIME 00:10:00

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/neuroscan32/Subject1_MP.cnt');

cfg.dataset = filename;
cfg.demean  = 'yes';
dataall     = ft_preprocessing(cfg);

cfg.trl     = [1 697 0;15335 16031 0;16032 16728 0];
data        = ft_preprocessing(cfg);

% the following assertion is needed due to numerical issues, the data are
% the same.
assert(all(abs(reshape(ft_preproc_baselinecorrect(dataall.trial{1}(:,1:697))'       - data.trial{1}',[],1))<1e-11));
assert(all(abs(reshape(ft_preproc_baselinecorrect(dataall.trial{1}(:,15335:16031))' - data.trial{2}',[],1))<1e-11));
assert(all(abs(reshape(ft_preproc_baselinecorrect(dataall.trial{1}(:,16032:16728))' - data.trial{3}',[],1))<1e-11));

assert(~isequal(data.trial{1},data.trial{3}));

