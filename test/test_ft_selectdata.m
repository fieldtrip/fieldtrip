function test_ft_selectdata

% TEST test_ft_selectdata
% TEST ft_selectdata ft_selectdata_old ft_selectdata_new ft_appendfreq

clear freq*

% make some dummy frequency structures
freq1.label = {'1' '2'};
freq1.freq  = 1:10;
freq1.time  = 1:5;
freq1.dimord = 'chan_freq_time';
freq1.powspctrm = randn(2,10,5);

cfg = [];
cfg.parameter = 'powspctrm';
freq2  = ft_appendfreq(cfg, freq1, freq1);
freq2  = rmfield(freq2, 'cfg');
freq2a = ft_selectdata(freq1, freq1, 'param', 'powspctrm'); % this should append the power spectrum
assert(isequal(freq2, freq2a));

freq4a = ft_selectdata(freq2, freq2, 'param', 'powspctrm');
assert(isequal(size(freq4a.powspctrm), [4 2 10 5]));

clear freq*

freq3.label = {'1' '2'};
freq3.freq  = 1:10;
freq3.dimord = 'chan_freq';
freq3.powspctrm = randn(2,10);

cfg = [];
cfg.parameter = 'powspctrm';
freq4  = ft_appendfreq(cfg, freq3, freq3);
freq4  = rmfield(freq4, 'cfg');
freq4a = ft_selectdata(freq3, freq3, 'param', 'powspctrm');  % this should append the power spectrum
assert(isequal(freq4, freq4a));

timelock1 = [];
timelock1.label = {'1' '2'};
timelock1.time  = 1:5;
timelock1.dimord = 'chan_time';
timelock1.avg = randn(2,5);

cfg = [];
cfg.channel = 1;
timelock1a = ft_selectdata(cfg, timelock1);
assert(isequal(size(timelock1a.avg), [1 5]));

cfg = [];
timelock2 = ft_appendtimelock(cfg, timelock1, timelock1, timelock1);

cfg = [];
cfg.channel = 1;
timelock2a = ft_selectdata(cfg, timelock2);
assert(isequal(size(timelock2a.trial), [3 1 5]));

cfg = [];
cfg.trials = [1 2];
timelock2b = ft_selectdata(cfg, timelock2);
assert(isequal(size(timelock2b.trial), [2 2 5]));

% The one that follows is a degenerate case. By selecting only one trial,
% the output is not really trial-based any more, but still contains one trial.
cfg = [];
cfg.trials = 1;
timelock2c = ft_selectdata(cfg, timelock2);
assert(isequal(size(timelock2c.trial), [1 2 5]));
% assert(isequal(size(timelock2c.trial), [2 5]));
