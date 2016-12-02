function test_bug3214

% WALLTIME 00:20:00
% MEM 2gb

%% try to simulate some data

data20 = [];
data20.label = {'1', '2', '3'};
for i=1:20
  data20.trial{i} = randn(3,2000);
  data20.time{i} = (1:2000)/1000;
end

data25 = [];
data25.label = {'1', '2', '3'};
for i=1:25
  data25.trial{i} = randn(3,2000);
  data25.time{i} = (1:2000)/1000;
end

%% this is what was done according to the bug report

cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.keeptrials = 'yes';
cfg.foi          = 0:2:40;  %% 20 freqs!
cfg.t_ftimwin    = ones(length(cfg.foi),1).*.5;
cfg.toi          = 0:0.05:2;
tfr20 = ft_freqanalysis(cfg, data20); %% data has 20 trials!
tfr25 = ft_freqanalysis(cfg, data25); %% data has 25 trials!

%% and here is the bug

cfg = [];
cfg.trials = 1:5;
sel20 = ft_selectdata(cfg, tfr20); % it was here
sel25 = ft_selectdata(cfg, tfr25); % but not here

%% this is to check another condition

cfg = [];
cfg.output       = 'fourier';
cfg.channel      = 'all';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.keeptrials = 'yes';
freq20 = ft_freqanalysis(cfg, data20);
freq25 = ft_freqanalysis(cfg, data25);

cfg = [];
cfg.trials = 1:5;
sel20a = ft_selectdata(cfg, freq20);
sel25a = ft_selectdata(cfg, freq25);

% try to make it more confusing by transposing the vector
freq20.cumsumcnt = freq20.cumsumcnt';
freq20.cumtapcnt = freq20.cumtapcnt';

cfg = [];
cfg.trials = 1:5;
sel20b = ft_selectdata(cfg, freq20);
sel25b = ft_selectdata(cfg, freq25);
