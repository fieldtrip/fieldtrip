function test_ft_timelockgrandaverage

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_timelockgrandaverage ft_timelockanalysis

% This functions tests the new implementation of ft_timelockgrandaverage. The new
% functionality includes the use of a cfg.parameter, and allows for 'chan' data to be
% averaged/combined within this function. The first part just tests whether it runs
% through for the most commonly used applications. It doesn't test the toilim
% selection.

%%
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

cfg.parameter = 'avg';
grandavg2 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.parameter = 'stat';
grandavg3 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.keepindividual = 'yes';
grandavg4 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

%%
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

cfg.parameter = 'avg';
grandavg2 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.parameter = 'stat';
grandavg3 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

cfg.keepindividual = 'yes';
grandavg4 = ft_timelockgrandaverage(cfg, timelock1, timelock2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the following code checks the functionality of the cfg.method = 'within'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = [];
data.trial = cat(2,repmat({[1 1]},[1 4]),repmat({[2 2]},[1 3]),repmat({[3 3]},[1 2]));
data.time  = repmat({[1 2]},[1 9]);
data.label = {'chan01'};

cfg = [];
cfg.normalizevar = 'N-1';
timelock = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.normalizevar = 'N-1';
cfg.trials = 1:4;
tlck1 = ft_timelockanalysis(cfg, data);
cfg.trials = 5:7;
tlck2 = ft_timelockanalysis(cfg, data);
cfg.trials = 8:9;
tlck3 = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.method = 'within';
cfg.normalizevar = 'N-1';
grandavg = ft_timelockgrandaverage(cfg, tlck1, tlck2, tlck3);

timelock = rmfield(timelock, 'cfg');
grandavg = rmfield(grandavg, 'cfg');

assert(isalmostequal(timelock, grandavg, 'reltol', 1000*eps));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the following code checks the functionality of data with nans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create some data containing nans, these would occur if there are artifacts
clear data_sub1 data_sub2 data_sub3

% one trial with nans
data_sub1.label = {'chan1', 'chan2'};
data_sub1.time{1} = 1:5;
data_sub1.fsample = 1;
data_sub1.trial{1} = rand(2,5);
data_sub1.trial{1}(1,1) = nan;
data_sub1.trialinfo = [1];

% two trials, one has nans
data_sub2.label = {'chan1'; 'chan2'};
data_sub2.time{1} = 1:5;
data_sub2.time{2} = 1:5;
data_sub2.fsample = 1;
data_sub2.trial{1} = rand(2,5);
data_sub2.trial{2} = rand(2,5);
data_sub2.trial{2}(1,1:2) = nan;
data_sub2.trialinfo = [1;1];

% three trials, one has nans
data_sub3.label = {'chan1'; 'chan2'};
data_sub3.time{1} = 1:5;
data_sub3.time{2} = 1:5;
data_sub3.time{3} = 1:5;
data_sub3.fsample = 1;
data_sub3.trial{1} = rand(2,5);
data_sub3.trial{2} = rand(2,5);
data_sub3.trial{3} = rand(2,5);
data_sub3.trial{3}(1,1:3) = nan;
data_sub3.trialinfo = [1;1;1];

% also append all trials, this allows an alternative to computing within-subject variance
cfg = [];
cfg.keeptrials = 'yes';
timelock_keep1 = ft_timelockanalysis(cfg, data_sub1);
timelock_keep2 = ft_timelockanalysis(cfg, data_sub2);
timelock_keep3 = ft_timelockanalysis(cfg, data_sub3);

cfg = [];
cfg.parameter = 'trial';
timelock_keep12  = ft_appendtimelock(cfg, timelock_keep1, timelock_keep2);
timelock_keep123 = ft_appendtimelock(cfg, timelock_keep1, timelock_keep2, timelock_keep3);

assert(isequal(size(timelock_keep12.trial), [3 2 5]));
assert(isequal(size(timelock_keep123.trial), [6 2 5]));

%%
% timelock

cfg = [];
% cfg.nanmean = 'yes';        % this is the default
% cfg.normalizevar = 'N-1';   % this is the default
timelock_sub1 = ft_timelockanalysis(cfg, data_sub1);
timelock_sub2 = ft_timelockanalysis(cfg, data_sub2);
timelock_sub3 = ft_timelockanalysis(cfg, data_sub3);

% check
assert( any(isnan(timelock_sub1.avg(:)))); % since this is based on a single trial with nans, there are still nans
assert(~any(isnan(timelock_sub2.avg(:)))); % no nan in the average
assert(~any(isnan(timelock_sub3.avg(:)))); % no nan in the average

assert( any(isnan(timelock_sub1.var(:)))); % since this is based on a single trial with nans, there are still nans
assert( any(isnan(timelock_sub2.var(:)))); % since this is based on two trials where one had nans, there are still nans
assert(~any(isnan(timelock_sub3.var(:)))); % no nan in the variance

% ensure that the nans are at the expected place
assert(isequal(isnan(timelock_sub1.avg), [1 0 0 0 0; 0 0 0 0 0]));
assert(isequal(isnan(timelock_sub1.var), [1 1 1 1 1; 1 1 1 1 1])); % since only a single trial

% ensure that degrees of freedom are correct, these are sum(~isnan(..)) of the concatenated trials
assert(isequal(timelock_sub1.dof, [0 1 1 1 1; 1 1 1 1 1]));
assert(isequal(timelock_sub2.dof, [1 1 2 2 2; 2 2 2 2 2]));
assert(isequal(timelock_sub3.dof, [2 2 2 3 3; 3 3 3 3 3]));

%%
% grand average

cfg = [];
cfg.keepindividual = 'yes';
grandavg_keep12  = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2);
grandavg_keep123 = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2, timelock_sub3);

assert(isequal(size(grandavg_keep12.individual), [2 2 5]));
assert(isequal(size(grandavg_keep123.individual), [3 2 5]));

%%

cfg = [];
cfg.keepindividual = 'no';
cfg.nanmean = 'yes';        % this is the default
cfg.normalizevar = 'N-1';   % this is the default

cfg.method = 'across';
grandavg_across12  = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2);
grandavg_across123 = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2, timelock_sub3);

cfg.method = 'within';
grandavg_within12  = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2);
grandavg_within123 = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2, timelock_sub3);

% check
assert(~any(isnan(grandavg_across12.avg(:))));
assert( any(isnan(grandavg_across12.var(:))));

assert(~any(isnan(grandavg_across123.avg(:))));
assert(~any(isnan(grandavg_across123.var(:))));

assert(~any(isnan(grandavg_within12.avg(:))));
assert( any(isnan(grandavg_within12.var(:)))); % since subject 1 variance is all nan

assert(isequal(grandavg_within12.dof, squeeze(sum(isfinite(timelock_keep12.trial),1))));
assert(isalmostequal(grandavg_within12.avg, squeeze(nanmean(timelock_keep12.trial,1)), 'reltol', 100*eps));

assert(~any(isnan(grandavg_within123.avg(:))));
assert(~any(isnan(grandavg_within123.var(:)))); % since subject 1 variance is all nan

assert(isequal(grandavg_within123.dof, squeeze(sum(isfinite(timelock_keep123.trial),1))));
assert(isalmostequal(grandavg_within123.avg, squeeze(nanmean(timelock_keep123.trial,1)), 'reltol', 100*eps));

%%

cfg = [];
cfg.keepindividual = 'no';
cfg.nanmean = 'no';         % DIFFERENCE TO ABOVE
cfg.normalizevar = 'N-1';   % this is the default

cfg.method = 'across';
grandavg_across12  = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2);
grandavg_across123 = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2, timelock_sub3);

cfg.method = 'within';
grandavg_within12  = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2);
grandavg_within123 = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2, timelock_sub3);

% check
assert(any(isnan(grandavg_across12.avg(:))));
assert(any(isnan(grandavg_across12.var(:))));

assert(any(isnan(grandavg_across123.avg(:))));
assert(any(isnan(grandavg_across123.var(:))));

assert(any(isnan(grandavg_within12.avg(:))));
assert(any(isnan(grandavg_within12.var(:))));

assert(isequal(grandavg_within12.dof, squeeze(sum(isfinite(timelock_keep12.trial),1))));

assert(any(isnan(grandavg_within123.avg(:))));
assert(any(isnan(grandavg_within123.var(:))));

assert(isequal(grandavg_within123.dof, squeeze(sum(isfinite(timelock_keep123.trial),1))));