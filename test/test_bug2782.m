function test_bug2782

% WALLTIME 00:10:00
% MEM 1gb


%% this is the initial problem

timelock = [];
timelock.label = {'1', '2', '3'};
timelock.time = (1:1000)/1000;
timelock.avg = randn(3,1000);
timelock.trial = randn(1,3,1000);
timelock.dimord = 'rpt_chan_time';
timelock.trialinfo = 1:11;

cfg = [];
cfg.trials = 1;
tmp = ft_selectdata(cfg,timelock);

% this should not change
assert(size(tmp.trialinfo,1)==1);
assert(size(tmp.trialinfo,2)==11);

%% this is according to comment 10

raw = [];
raw.label = {'1', '2', '3'};
raw.time{1}  = (1:1000)/1000;
raw.trial{1} = randn(3,1000);
raw.sampleinfo = [1 1000];


cfg = [];
cfg.trials = 1;
tmp = ft_selectdata(cfg,raw);

% this should not change
assert(size(tmp.sampleinfo,1)==1);
assert(size(tmp.sampleinfo,2)==2);

%% this is according to comment 7

raw = [];
raw.label = {'1', '2', '3'};
for i=1:5
  raw.time{i}  = (1:1000)/1000;
  raw.trial{i} = i + 0.1*randn(3,1000);
  raw.sampleinfo(i,:) = [1 1000] + (i-1)*1000;
end

cfg = [];
cfg.trials = 5;
cfg.latency = [0.1 0.3];
tmp = ft_selectdata(cfg, raw);

assert(tmp.time{1}(  1)==0.1);
assert(tmp.time{1}(end)==0.3);
assert(abs(mean(tmp.trial{1}(:))-5)<0.1);


