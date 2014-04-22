function test_bug2508

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_bug2508 ft_selectdata ft_selectdata_new

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check the averaging over dimensions for a simple timelock structure

timelock = [];
timelock.label = {'1', '2', '3'};
timelock.time = 1:5;
timelock.avg = randn(3,5);
timelock.dimord = 'chan_time';

cfg = [];
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
output = ft_selectdata(cfg, timelock);

assert(length(output.label)==3, 'error in output labels');
assert(length(output.time)==5,  'error in output time');

cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'no';
output = ft_selectdata(cfg, timelock);

assert(length(output.label)==1, 'error in output labels');
assert(length(output.time)==5,  'error in output time');

cfg = [];
cfg.avgoverchan = 'no';
cfg.avgovertime = 'yes';
output = ft_selectdata(cfg, timelock);

assert(length(output.label)==3, 'error in output labels');
assert(length(output.time)==1,  'error in output time');

cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
output = ft_selectdata(cfg, timelock);

assert(length(output.label)==1, 'error in output labels');
assert(length(output.time)==1,  'error in output time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% continue with timelock, now with repetitions

timelock = [];
timelock.label = {'1', '2', '3'};
timelock.time = 1:5;
timelock.trial = randn(2,3,5);
timelock.dimord = 'rpt_chan_time';

cfg = [];
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverrpt = 'no';
output = ft_selectdata(cfg, timelock);

assert(isequal(size(output.trial), [2 3 5]), 'error in output data size');
assert(length(output.label)==3, 'error in output labels');
assert(length(output.time)==5,  'error in output time');

cfg = [];
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverrpt = 'yes';
output = ft_selectdata(cfg, timelock);

assert(isequal(size(output.trial), [3 5]), 'error in output data size'); % note the confusing average field name output.trial
assert(length(output.label)==3, 'error in output labels');
assert(length(output.time)==5,  'error in output time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now see what happens when some fields are not kept

cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverrpt  = 'yes';
cfg.keeprptdim  = 'yes';
cfg.keeptimedim = 'yes';
cfg.keepchandim = 'yes';
output = ft_selectdata(cfg, timelock);

assert(strcmp(output.dimord, 'rpt_chan_time'), 'incorrect dimord');
assert(numel(output.trial)==1, 'error in output data size'); % note the confusing average field name output.trial
assert(length(output.label)==1, 'error in output labels');
assert(length(output.time)==1,  'error in output time');

cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverrpt  = 'yes';
cfg.keeprptdim  = 'no';
cfg.keeptimedim = 'yes';
cfg.keepchandim = 'yes';
output = ft_selectdata(cfg, timelock);

assert(strcmp(output.dimord, 'chan_time'), 'incorrect dimord');
assert(numel(output.trial)==1, 'error in output data size'); % note the confusing average field name output.trial
assert(length(output.label)==1, 'error in output labels');
assert(length(output.time)==1,  'error in output time');

cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverrpt  = 'yes';
cfg.keeprptdim  = 'no';
cfg.keeptimedim = 'no';
cfg.keepchandim = 'yes';
output = ft_selectdata(cfg, timelock);

assert(strcmp(output.dimord, 'chan'), 'incorrect dimord');
assert(numel(output.trial)==1, 'error in output data size'); % note the confusing average field name output.trial
assert(~isfield(output, 'time'),  'output should not have time field');
assert( isfield(output, 'label'), 'output should have time field');

cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverrpt  = 'yes';
cfg.keeprptdim  = 'no';
cfg.keeptimedim = 'yes';
cfg.keepchandim = 'no';
output = ft_selectdata(cfg, timelock);

assert(strcmp(output.dimord, 'time'), 'incorrect dimord');
assert(numel(output.trial)==1, 'error in output data size'); % note the confusing average field name output.trial
assert( isfield(output, 'time'),  'output should have time field');
assert(~isfield(output, 'label'), 'output should not have time field');
