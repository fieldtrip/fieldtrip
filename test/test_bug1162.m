function test_bug1162

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_postamble ft_postamble_history ft_freqgrandaverage

timelock1 = [];
timelock1.label = {'1' '2'};
timelock1.time  = 1:5;
timelock1.dimord = 'rpt_chan_time';
timelock1.avg = randn(1,2,5);
timelock1.cfg = 'this is number 1';

timelock2 = [];
timelock2.label = {'1' '2'};
timelock2.time  = 1:5;
timelock2.dimord = 'rpt_chan_time';
timelock2.avg = randn(1,2,5);
timelock2.cfg = 'this is number 2';

cfg = [];
grandavg = ft_timelockgrandaverage(cfg, timelock1, timelock2);

% grandavg.cfg.previous{1} should be 'this is number 1'
% grandavg.cfg.previous{2} should be 'this is number 1'
disp(grandavg.cfg.previous);

assert(iscell(grandavg.cfg.previous));
assert(~isempty(grandavg.cfg.previous));
assert(~isempty(grandavg.cfg.previous{1}));
assert(~isempty(grandavg.cfg.previous{2}));

