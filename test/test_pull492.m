function test_pull492

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

%% 

nchan = 10;
ntime = 1000;
ntrial = 42;

data = [];
for i=1:nchan
  data.label{i} = num2str(i);
end

for i=1:ntrial
  data.time{i} = (1:ntime)/ntime;
  data.trial{i} = randn(nchan, ntime);
  data.trial{i}(1,:) = i;  % replace values with the trial number
end


%% 


cfg = [];
cfg.trials = setdiff(1:42, 3);
cfg.keeptrials = 'yes';
timelock3 = ft_timelockanalysis(cfg, data);
dum = timelock3.trial(:,1,1);
assert(isequal(dum', cfg.trials));

cfg = [];
cfg.trials = setdiff(1:42, [3, 6, 9]);
cfg.keeptrials = 'yes';
timelock3 = ft_timelockanalysis(cfg, data);
dum = timelock3.trial(:,1,1);
assert(isequal(dum', cfg.trials));

cfg = [];
cfg.trials = 'all';
cfg.keeptrials = 'yes';
timelock = ft_timelockanalysis(cfg, data);

dum = timelock.trial(:,1,1);
assert(isequal(dum', 1:42));

%%

cfg = [];
cfg.method = 'wavelet';
cfg.keeptrials = 'yes';
cfg.toi = 0:0.02:1;
cfg.foi = 1:100; 
freq = ft_freqanalysis(cfg, data);
% set the first channel-time-frequency value equal to the trial number
freq.powspctrm(:,1,1,1) = 1:ntrial;

%%

layout = [];
for i=1
  % only the first channel is sufficient
  layout.label{i}   = num2str(i);
  layout.width(i)   = 1;
  layout.height(i)  = 2;
  layout.pos(i,:)   = [0.5 * (ntime+1)/ntime, 0];
end
% ensure that it does not autoscale in the default spherical head
layout.outline = {};
layout.mask    = {};

%%

cfg = [];
cfg.layout = layout;
cfg.trials = [1 4 7];
ft_multiplotER(cfg, timelock);

axis on
axis([0 1 -45 45])

% the last one corresponds to channel 1
c = get(gca, 'Children');
ydata = get(c(end), 'Ydata');

assert(ydata(1) == mean(cfg.trials));

%%

cfg = [];
cfg.layout = layout;
cfg.trials = [1 4 7];
ft_multiplotTFR(cfg, freq);

% the last one corresponds to channel 1
c = get(gca, 'Children');
cdata = get(c(end), 'CData');
assert(cdata(1,1) == mean(cfg.trials));

