function test_timelockNaN

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_timelock ft_timelockgrandaverage

% This function tests whether averages and variances are calculated
% correctly when using ft_timelock and ft_timelockgrandaverage with data
% containing nans.

%%
% create some data containing nans, these would occur if there are artifacts
clear data_sub1 data_sub2 data_sub3

% one trial with nans
data_sub1.label = {'chan1', 'chan2'};
data_sub1.time{1} = 1:5;
data_sub1.fsample = 1;
data_sub1.trial{1} = rand(2,5);
data_sub1.trial{1}(1,1:2) = nan;
data_sub1.trialinfo = [1];

% two trials, one has nans
data_sub2.label = {'chan1'; 'chan2'};
data_sub2.time{1} = 1:5;
data_sub2.time{2} = 1:5;
data_sub2.fsample = 1;
data_sub2.trial{1} = rand(2,5);
data_sub2.trial{2} = rand(2,5);
data_sub2.trial{2}(1,2:3) = nan;
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
data_sub3.trial{3}(1,3:4) = nan;
data_sub3.trialinfo = [1;1;1];

%%
% timelock
cfg = [];
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
assert(isequal(isnan(timelock_sub1.avg), [1 1 0 0 0; 0 0 0 0 0]));
assert(isequal(isnan(timelock_sub1.var), [1 1 1 1 1; 1 1 1 1 1])); % since only a single trial
assert(isequal(isnan(timelock_sub2.var), [0 1 1 0 0; 0 0 0 0 0]));

% ensure that degrees of freedom are correct, these are sum(~isnan(..)) of the concatenated trials
assert(isequal(timelock_sub1.dof, [0 0 1 1 1; 1 1 1 1 1]));
assert(isequal(timelock_sub2.dof, [2 1 1 2 2; 2 2 2 2 2]));
assert(isequal(timelock_sub3.dof, [3 3 2 2 3; 3 3 3 3 3]));

%%
% grand average
cfg = [];
cfg.method = 'across';
grandavg_across12  = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2);
grandavg_across123 = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2, timelock_sub3);
cfg.method = 'within';
grandavg_within12  = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2);
grandavg_within123 = ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2, timelock_sub3);

% check
assert(~any(isnan(grandavg_across12.avg(:))));
assert(~any(isnan(grandavg_across123.avg(:))));
assert( any(isnan(grandavg_across12.var(:))));  % since partially nan in one of the two subjects
assert(~any(isnan(grandavg_across123.var(:))));

% FIXME these do not seem correct, but I am also not 100% sure what to expect
% assert(~any(isnan(grandavg_within12.avg(:))));
% assert(~any(isnan(grandavg_within123.avg(:))));
% assert( any(isnan(grandavg_within12.var(:))));
% assert(~any(isnan(grandavg_within123.var(:))));
