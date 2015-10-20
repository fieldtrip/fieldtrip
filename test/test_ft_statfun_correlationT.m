function test_ft_statfun_correlationT

% simulate simple single-subject timelock structures
data_brain = [];
data_brain.trial = [1:10]'; % increasing 
data_brain.dimord = 'rpt_chan_time';
data_brain.time = 1;
data_brain.label = {'1'};

data_behav = data_brain;
data_behav.trial = data_brain.trial*-1000+50; % add scaling difference

% compute statistics with correlationT
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_correlationT';
cfg.numrandomization = 100;

n1 = 10;    % n1 is the number of trials
design              = zeros(2, n1 * 2);
design(1,1:n1)      = 1;
design(1,(n1 + 1):(n1 * 2)) = 2;
design(2, :)        = [1:n1 1:n1];
cfg.design           = design;

cfg.ivar             = 1;
cfg.uvar             = 2;
stat = ft_timelockstatistics(cfg, data_brain, data_behav);

assert(isequal(stat.rho, -1));


% simulate simple multiple subjects timelock structures
data_brain = [];
data_behav = [];
for j=1:10
data_brain{j}.avg = j; % increasing 
data_brain{j}.dimord = 'chan_time';
data_brain{j}.time = 1;
data_brain{j}.label = {'1'};

data_behav{j} = data_brain{j};
data_behav{j}.avg = data_brain{j}.avg*-1000+50; % add scaling difference
end

% compute statistics with correlationT
cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_correlationT';
cfg.numrandomization = 100;

n1 = 10;    % n1 is the number of subjects
design              = zeros(2, n1 * 2);
design(1,1:n1)      = 1;
design(1,(n1 + 1):(n1 * 2)) = 2;
design(2, :)        = [1:n1 1:n1];
cfg.design           = design;

cfg.ivar             = 1;
cfg.uvar             = 2;
stat = ft_timelockstatistics(cfg, data_brain{:}, data_behav{:});

assert(isequal(stat.rho, -1));