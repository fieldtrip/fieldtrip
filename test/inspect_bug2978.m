function inspect_bug2978

ntrial = 5;

%%

data = [];
data.label = {'1', '2', '3'};
for i=1:ntrial
  t = (1:1000)/1000;
  s = sin(i*2*pi*t); % make a sine with 1, 2, 3, ... Hz
  data.time{i} = t;
  data.trial{i} = 0.1*randn(3,1000);
  data.trial{i}(1,:) = s;  % insert it in the first channel
end

%%

cfg = [];
timelock = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.layout = 'vertical';
figure
ft_multiplotER(cfg, timelock); % it should show a mix


cfg = [];
cfg.trials = 1;
timelock1 = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.trials = 2;
timelock2 = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.layout = 'vertical';
figure
ft_multiplotER(cfg, timelock1, timelock2); % it should show 1 and 2 Hz

%%

cfg = [];
cfg.keeptrials = 'yes';
timelock_all = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.avgoverrpt = 'yes';
cfg.trials = 1;
timelock1 = ft_selectdata(cfg, timelock_all);

figure
plot(timelock1.time, timelock1.avg) % this stays the same, i.e. a frequency mix
figure
plot(timelock1.time, timelock1.trial) % this should show 1 Hz




%%
% THESE TWO ARE STILL FAILING ON 19 May 2017

if false
  
  cfg = [];
  cfg.layout = 'vertical';
  cfg.trials = 1;
  figure
  ft_multiplotER(cfg, data); % it should show 1Hz
  
  cfg.trials = 2;
  figure
  ft_multiplotER(cfg, data); % it should show 2Hz
  
end

%%

cfg = [];
cfg.layout = 'vertical';
cfg.method = 'summary';
ft_rejectvisual(cfg, data);
