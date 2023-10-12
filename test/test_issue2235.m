function test_issue2235

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY data2bids
% DATA no

%% make an EEG dataset in memory

fsample = 1000;

data = [];
data.label = {'1', '2', '3'};
data.time{1} = ((1:10000)-1)/fsample;
data.trial{1} = randn(3, 10000);

%% write it to bids, with a single event

onset = 0;
duration = 10;
type = {'trial'};

cfg = [];
cfg.bidsroot = tempname;
cfg.sub = '01';
cfg.ses = '01';
cfg.datatype = 'eeg';
cfg.task = 'eyesopen';
cfg.events = table(onset, duration, type);

data2bids(cfg, data)

%% clean up

disp('cleaning up')
rmdir(cfg.bidsroot, 's');
