function test_ft_statfun_correlationT

% A. Stolk, oct 2015

% simulate simple single-subject timelock structures
data_brain = [];
data_brain.trial = rand(10,1); % increasing
data_brain.dimord = 'rpt_chan_time';
data_brain.time = 1;
data_brain.label = {'1'};

data_behav = data_brain;
data_behav.trial = rand(10,1)+1000; % add offset (which should be not allowed to bias the results)

% compute statistics with correlationT
cfg = [];
cfg.statistic        = 'ft_statfun_correlationT';
cfg.method           = 'montecarlo';
cfg.resampling       = 'bootstrap'; % this parameter is crucial, the default (permutation) is prone to systematic across-condition bias
cfg.numrandomization = 1000;

n1 = 10;    % n1 is the number of trials
design              = zeros(2, n1 * 2);
design(1,1:n1)      = 1;
design(1,(n1 + 1):(n1 * 2)) = 2;
design(2, :)        = [1:n1 1:n1];
cfg.design           = design;

cfg.ivar             = 1;
cfg.uvar             = 2;
stat1 = ft_timelockstatistics(cfg, data_brain, data_behav);

% simulate simple multiple subjects timelock structures
data_brain = [];
data_behav = [];
for j=1:10
  data_brain{j}.avg = rand; % increasing
  data_brain{j}.dimord = 'chan_time';
  data_brain{j}.time = 1;
  data_brain{j}.label = {'1'};
  
  data_behav{j} = data_brain{j};
  data_behav{j}.avg = rand+1000; % add scaling difference
end

% compute statistics with correlationT
cfg = [];
cfg.statistic        = 'ft_statfun_correlationT';
cfg.method           = 'montecarlo';
cfg.resampling       = 'bootstrap'; % this parameter is crucial, the default (permutation) is prone to systematic across-condition bias
cfg.numrandomization = 1000;

n1 = 10;    % n1 is the number of subjects
design              = zeros(2, n1 * 2);
design(1,1:n1)      = 1;
design(1,(n1 + 1):(n1 * 2)) = 2;
design(2, :)        = [1:n1 1:n1];
cfg.design           = design;

cfg.ivar             = 1;
cfg.uvar             = 2;
stat2 = ft_timelockstatistics(cfg, data_brain{:}, data_behav{:});

assert(~stat1.prob<0.05 && ~stat2.prob<0.05); % probability-wise, this should fail once every millenium :)