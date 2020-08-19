function test_timelockNaN

% DEPENDENCY ft_timelock, ft_timelockgrandaverage

% This function tests whether averages and variances are calculated
% correctly when using ft_timelock and ft_timelockgrandaverage with data
% containing nans.

% create some data containing nans
clear data_sub1 data_sub2
data_sub1.label={'chan1'; 'chan2'};
data_sub1.time{1}=1:5;
data_sub1.time{2}=1:5;
data_sub1.fsample=1;
data_sub1.trial{1}=rand(2,5);
data_sub1.trial{2}=rand(2,5);
data_sub1.trial{2}(1,3:5)=nan;
data_sub1.trialinfo=[1;1];

data_sub2.label={'chan1', 'chan2'};
data_sub2.time{1}=1:5;
data_sub2.fsample=1;
data_sub2.trial{1}=rand(2,5);
data_sub2.trial{1}(1,1:2)=nan;

% timelock
cfg=[];
timelock_sub1=ft_timelockanalysis(cfg, data_sub1);
timelock_sub2=ft_timelockanalysis(cfg, data_sub2);
% check
assert(isempty(find(isnan(timelock_sub1.avg))));
assert(sum(isnan(timelock_sub1.var)==[0 0 1 1 1;0 0 0 0 0], 'all')==10);
assert(sum(isnan(timelock_sub2.avg)==[1 1 0 0 0;0 0 0 0 0], 'all')==10);
assert(sum(isnan(timelock_sub2.var)==ones(2,5), 'all')==10);

% grand average
cfg=[];
GA=ft_timelockgrandaverage(cfg, timelock_sub1, timelock_sub2);
% check
assert(isempty(find(isnan(GA.avg))));
assert(sum(isnan(GA.var)==[1 1 0 0 0;0 0 0 0 0], 'all')==10);