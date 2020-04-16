function test_ft_regressconfound

% MEM 2gb
% WALLTIME 00:10:00

% DEPENDENCY ft_regressconfound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq1 = [];
freq1.label = {'1' '2'};
freq1.freq  = 1:10;
freq1.dimord = 'rpt_chan_freq';
freq1.powspctrm = randn(20,2,10);

cfg = [];
cfg.confound = randn(20,3);
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'}; % this needs stat toolbox
freq1_out = ft_regressconfound(cfg, freq1);
assert(isfield(freq1_out, 'powspctrm'));
assert(~isfield(freq1_out, 'beta'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq2 = [];
freq2.label = {'1' '2'};
freq2.freq  = 1:10;
freq2.time  = 1:5;
freq2.dimord = 'rpt_chan_freq_time';
freq2.powspctrm = randn(20,2,10,5);

cfg = [];
cfg.confound = randn(20,3);
cfg.output = 'beta';
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
freq2_out = ft_regressconfound(cfg, freq2);
assert(isfield(freq2_out, 'beta'));
assert(~isfield(freq2_out, 'powspctrm'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock = [];
timelock.label = {'1' '2'};
timelock.time  = 1:5;
timelock.dimord = 'rpt_chan_time';
timelock.trial = randn(20,2,5);
timelock.avg = randn(2,5);

cfg = [];
cfg.confound = randn(20,3);
cfg.output = 'model';
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
timelock_out = ft_regressconfound(cfg, timelock);
assert(isfield(timelock_out, 'model'));
assert(~isfield(timelock_out, 'trial'));

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
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
source_out = ft_regressconfound(cfg, source);
assert(isfield(source_out, 'pow'));
assert(~isfield(source_out, 'beta'));

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
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
timelock2_out = ft_regressconfound(cfg, timelock2);
assert(isfield(timelock2_out, 'trial'));
assert(~isfield(timelock2_out, 'beta'));

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
cfg.reject = [1:3];
%cfg.Ftest = {'1','2','3'};
source2_out = ft_regressconfound(cfg, source2);
assert(isfield(source2_out, 'pow'));
assert(~isfield(source2_out, 'beta'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source3 = [];
source3.pos = randn(100,3);
source3.inside = true(100,1);
source3.pow = randn(100,25);
source3.dimord = 'pos_rpt';

cfg = [];
cfg.confound = randn(25,2);
cfg.normalize = 'no';
source3_out = ft_regressconfound(cfg, source3);

betas = cfg.confound \ transpose(source3.pow);
desired = transpose(transpose(source3.pow) - cfg.confound * betas);

assert(isequal(source3_out.pow, desired));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source4 = [];
source4.pos = randn(100,3);
source4.inside = true(100,1);
source4.pow = randn(100,25,10);
source4.time = 1:10;
source4.dimord = 'pos_rpt_time';

cfg = [];
cfg.confound = randn(25,2);
cfg.normalize = 'no';
source4_out = ft_regressconfound(cfg, source4);

assert(isequal(size(source4_out.pow), size(source4.pow)));
for t = 1:10
  betas = cfg.confound \ transpose(source4.pow(:,:,t));
  desired = transpose(transpose(source4.pow(:,:,t)) - cfg.confound * betas);
  assert(isequal(source4_out.pow(:,:,t), desired));
end


