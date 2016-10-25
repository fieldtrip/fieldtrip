function test_bug2965

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_bug2965
% TEST ft_checkdata

%%

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.pad = 'maxperlen';
cfg.feedback = 'no';

timelock = [];
timelock.time = (1:1000)/1000;
timelock.avg = randn(1,64,1000);
timelock.sampleinfo = [1 1000];
timelock.trialinfo = [1 2 3];
timelock.dimord = 'rpt_chan_time';
for i=1:64
  timelock.label{i} = num2str(i);
end

freq = ft_freqanalysis(cfg, timelock);
assert(isequal(size(freq.fourierspctrm), [1 64 501]))
assert(isfield(freq, 'trialinfo'));

%%
% other fieldname with dimord=rpt_chan_time

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.pad = 'maxperlen';
cfg.feedback = 'no';

timelock = [];
timelock.time = (1:1000)/1000;
timelock.whatever = randn(1,64,1000);
timelock.sampleinfo = [1 1000];
timelock.trialinfo = [1 2 3];
timelock.dimord = 'rpt_chan_time';
for i=1:64
  timelock.label{i} = num2str(i);
end

freq = ft_freqanalysis(cfg, timelock);
assert(isequal(size(freq.fourierspctrm), [1 64 501]))
assert(isfield(freq, 'trialinfo'));

%%
% idem, different dimord

timelock.dimord = 'subj_chan_time';
freq = ft_freqanalysis(cfg, timelock);
assert(isequal(size(freq.fourierspctrm), [1 64 501]))
assert(isfield(freq, 'trialinfo'));

%%
% simple avg

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.pad = 'maxperlen';
cfg.feedback = 'no';

timelock = [];
timelock.time = (1:1000)/1000;
timelock.avg = randn(64,1000);
timelock.dimord = 'chan_time';
for i=1:64
  timelock.label{i} = num2str(i);
end

freq = ft_freqanalysis(cfg, timelock);
assert(isequal(size(freq.fourierspctrm), [1 64 501]))
assert(~isfield(freq, 'trialinfo'));

%%
% other fieldname with dimord=chan_time

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.pad = 'maxperlen';
cfg.feedback = 'no';

timelock = [];
timelock.time = (1:1000)/1000;
timelock.whatever = randn(64,1000);
timelock.dimord = 'chan_time';
for i=1:64
  timelock.label{i} = num2str(i);
end

freq = ft_freqanalysis(cfg, timelock);
assert(isequal(size(freq.fourierspctrm), [1 64 501]))
assert(~isfield(freq, 'trialinfo'));
