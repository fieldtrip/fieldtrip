function test_ft_regressconfound

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_regressconfound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq1 = [];
freq1.label = {'1' '2'};
freq1.freq  = 1:10;
freq1.dimord = 'rpt_chan_freq';
freq1.powspctrm = randn(20,2,10);

cfg = [];
cfg.confound = randn(20,3);
cfg.model = 'yes';
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'}; % this needs stat toolbox
freq1_out = ft_regressconfound(cfg, freq1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq2 = [];
freq2.label = {'1' '2'};
freq2.freq  = 1:10;
freq2.time  = 1:5;
freq2.dimord = 'rpt_chan_freq_time';
freq2.powspctrm = randn(20,2,10,5);

cfg = [];
cfg.confound = randn(20,3);
cfg.model = 'yes';
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
freq2_out = ft_regressconfound(cfg, freq2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock = [];
timelock.label = {'1' '2'};
timelock.time  = 1:5;
timelock.dimord = 'rpt_chan_time';
timelock.trial = randn(20,2,5);
timelock.avg = randn(2,5);

cfg = [];
cfg.confound = randn(20,3);
cfg.model = 'yes';
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
timelock_out = ft_regressconfound(cfg, timelock);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source = [];
source.freq = 10;
source.pos = randn(5000,3);
source.inside = find(source.pos(:,1) > 0)';
source.outside = find(source.pos(:,1) <= 0)';
for j=1:20
  source.trial(1,j).pow = randn(size(source.pos(:,1)));
end

cfg = [];
cfg.confound = randn(20,3);
cfg.model = 'yes';
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
source_out = ft_regressconfound(cfg, source);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock2 = [];
timelock2.label = {'1','2'};
timelock2.time  = 1:5;
timelock2.dimord = 'rpt_chan_time';
timelock2.trial = randn(20,2,5);
timelock2.trial(1,1,1:5) = NaN;
timelock2.trial(3,1,3:5) = NaN;

cfg = [];
cfg.confound = randn(20,3);
cfg.model = 'yes';
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
timelock2_out = ft_regressconfound(cfg, timelock2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source2 = [];
source2.freq = 10;
source2.pos = randn(5000,3);
source2.inside = find(source2.pos(:,1) > 0)';
source2.outside = find(source2.pos(:,1) <= 0)';
for j=1:20
  source2.trial(1,j).pow = randn(size(source2.pos(:,1)));
end
source2.trial(1,10).pow = NaN(size(source2.pos(:,1)));

cfg = [];
cfg.confound = randn(20,3);
cfg.model = 'yes';
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
source2_out = ft_regressconfound(cfg, source2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
