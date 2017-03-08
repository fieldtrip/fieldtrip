function test_ft_timelockgrandaverage

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_timelockgrandaverage

% this functions tests the new implementation of ft_timelockgrandaverage. the
% new functionality includes the use of a cfg.parameter, and allows for
% 'chan' data to be averaged/combined within this function. it just tests
% whether it runs through for the most commonly used applications. it
% doesn't test the toilim selection

% create some data
timelock1.label = {'chan1';'chan2'};
timelock1.time  = 1:5;
timelock1.dimord = 'chan_time';
timelock1.avg  = randn(2,5);
timelock1.stat  = randn(2,5);

timelock2 = timelock1;
timelock2.avg  = randn(2,5);
timelock2.stat  = randn(2,5);

cfg = [];
grandavg1 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.parameter = 'stat';
grandavg2 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.parameter = {'avg' 'stat'};
grandavg3 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.keepindividual = 'yes';
grandavg4 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

% create some data with different triallengths
timelock1.label = {'chan1';'chan2'};
timelock1.time  = 1:5;
timelock1.dimord = 'chan_time';
timelock1.avg  = randn(2,5);
timelock1.stat  = randn(2,5);

timelock2 = timelock1;
timelock2.time  = 2:4;
timelock2.avg  = randn(2,3);
timelock2.stat  = randn(2,3);

cfg = [];
grandavg1 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.parameter = 'stat';
grandavg2 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.parameter = {'avg' 'stat'};
grandavg3 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.keepindividual = 'yes';
grandavg4 = ft_timelockgrandaverage(cfg, timelock1, timelock2);


% the following code checks the functionality of the cfg.method = 'within'
data = [];
data.trial = cat(2,repmat({[1 1]},[1 4]),repmat({[2 2]},[1 3]),repmat({[3 3]},[1 2]));
data.time  = repmat({[1 2]},[1 9]);
data.label = {'chan01'};
data.fsample = 1;

tlck = ft_timelockanalysis([], data);

cfg = [];
cfg.trials = 1:4;
tlck1 = ft_timelockanalysis(cfg, data);
cfg.trials = 5:7;
tlck2 = ft_timelockanalysis(cfg, data);
cfg.trials = 8:9;
tlck3 = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.method = 'within';
cfg.normalizevar = 'N-1';
tlckall = ft_timelockgrandaverage(cfg, tlck1, tlck2, tlck3);

tlck    = rmfield(tlck, 'cfg');
tlckall = rmfield(tlckall, 'cfg');

assert(isalmostequal(tlck, tlckall));

