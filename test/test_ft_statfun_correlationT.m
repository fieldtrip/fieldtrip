function test_ft_statfun_correlationT

% WALLTIME 00:20:00
% MEM 3gb

% A. Stolk, oct 2015

%% simulate simple single-subject timelock structures
n1 = 20;    % n1 is the number of trials
data_brain = [];
data_brain.trial = rand(n1,2,4);
data_brain.dimord = 'rpt_chan_time';
data_brain.label = {'1';'2'};
data_brain.time = [1:4];

data_behav = rand(n1,1); % inserted in the cfg.ivar-th row of cfg.design

% compute statistics with correlationT
cfg = [];
cfg.statistic        = 'ft_statfun_correlationT';
cfg.method           = 'montecarlo';
cfg.numrandomization = 1000;

design(1,1:n1)       = data_behav;
cfg.design           = design;
cfg.ivar             = 1;
stat1 = ft_timelockstatistics(cfg, data_brain);

%% simulate multiple single-subject timelock structures
n1 = 20;    % n1 is the number of subjects
data_brain = [];
for j=1:n1
  data_brain{j}.avg = rand(2,4);
  data_brain{j}.dimord = 'chan_time';
  data_brain{j}.label = {'1';'2'};
  data_brain{j}.time = [1:4];
end

data_behav = rand(n1,1); % inserted in the cfg.ivar-th row of cfg.design

% compute statistics with correlationT
cfg = [];
cfg.statistic        = 'ft_statfun_correlationT';
cfg.method           = 'montecarlo';
cfg.numrandomization = 1000;

design(1,1:n1)       = data_behav;
cfg.design           = design;
cfg.ivar             = 1;
stat2 = ft_timelockstatistics(cfg, data_brain{:});

%% simulate multiple single-subject timefreq structures
n1 = 20;    % n1 is the number of subjects
data_brain = [];
for j=1:n1
  data_brain{j}.powspctrm = rand(2,3,4);
  data_brain{j}.dimord = 'chan_freq_time';
  data_brain{j}.label = {'1';'2'};
  data_brain{j}.freq = [1:3];
  data_brain{j}.time = [1:4];
end

data_behav = rand(n1,1); % inserted in the cfg.ivar-th row of cfg.design

% compute statistics with correlationT
cfg = [];
cfg.statistic        = 'ft_statfun_correlationT';
cfg.method           = 'montecarlo';
cfg.numrandomization = 1000;

design(1,1:n1)       = data_behav;
cfg.design           = design;
cfg.ivar             = 1;
stat3 = ft_freqstatistics(cfg, data_brain{:});

assert(~stat1.prob(1)<0.05 && ~stat2.prob(1)<0.05 && ~stat3.prob(1)<0.05); % probability-wise, this should fail once every millenium :)
