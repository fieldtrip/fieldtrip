function test_pull1304

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_connectivityanalysis ft_connectivity_corr

%%

load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull1304.mat'));

%%

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1; % i.e. no real smoothing
cfg.foilim = [1 100];
[Data_freq] = ft_freqanalysis(cfg, Data);

%%

cfg = [];
cfg.method = 'plv';
cfg.partchannel = 'Stimuli';
Data_stat = ft_connectivityanalysis(cfg, Data_freq);

%%

assert(~any(isnan(Data_stat.plvspctrm(:))));

assert(length(Data_stat.label)==size(Data_stat.plvspctrm,1)); % nchan
assert(length(Data_stat.label)==size(Data_stat.plvspctrm,2)); % nchan

