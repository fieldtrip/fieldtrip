function test_bug1658

% WALLTIME 00:20:00
% MEM 4gb

% TEST ft_sourcedescriptives
% TEST ft_selectdata

% test whether ft_sourceescriptives supports trial selection

% load some source data
load(dccnpath(fullfile('/home/common/matlab/fieldtrip/data/test/latest/source/meg','source_grid_mtmfft_fourier_trl_DICS_keepall_rawtrial_ctf151.mat')));

% remove the spatial filters
source.trial = rmfield(source.trial, {'filter' 'filterdimord', 'label'});

cfg = [];
cfg.trials = [1 2 3 4];
source2 = ft_selectdata(cfg, source);

sel = find(source.inside, 1, 'first');
clear dat1 dat2;

for k = 1:numel(cfg.trials)
  dat1(k,:,:) = source.trial(cfg.trials(k)).csd{sel};
end
dat2 = source2.csd{sel};
assert(isequal(dat1,dat2));
clear dat1 dat2;

for k = 1:numel(cfg.trials)
  dat1(k) = source.trial(cfg.trials(k)).noise(sel);
end
dat2 = source2.noise(sel,:);
assert(isequal(dat1,dat2));
clear dat1 dat2

for k = 1:numel(cfg.trials)
  dat1(k) = source.trial(cfg.trials(k)).pow(sel);
end
dat2 = source2.pow(sel,:);
assert(isequal(dat1,dat2));

cfg = [];
cfg.trials = [1 2 3 4];
cfg.keeptrials = 'yes';
source3 = ft_sourcedescriptives(cfg, source);

sel = find(source.inside, 1, 'first');
clear dat1 dat2;

for k = 1:numel(cfg.trials)
  dat1(k,:,:) = source.trial(cfg.trials(k)).csd{sel};
  dat2(k,:,:) = source3.trial(cfg.trials(k)).csd{sel};
end
assert(isequal(dat1,dat2));
clear dat1 dat2;

for k = 1:numel(cfg.trials)
  dat1(k) = source.trial(cfg.trials(k)).noise(sel);
  dat2(k) = source3.trial(cfg.trials(k)).noise(sel);
end
assert(isequal(dat1,dat2));
clear dat1 dat2

for k = 1:numel(cfg.trials)
  dat1(k) = source.trial(cfg.trials(k)).pow(sel);
  dat2(k) = source3.trial(cfg.trials(k)).pow(sel);
end
assert(isequal(dat1,dat2));
