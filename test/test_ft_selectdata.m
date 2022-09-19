function test_ft_selectdata

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_selectdata ft_selectdata_old ft_selectdata_new ft_appendfreq

timelock1 = [];
timelock1.label = {'1' '2'};
timelock1.time  = 1:5;
timelock1.dimord = 'chan_time';
timelock1.avg = randn(2,5);

cfg = [];
cfg.channel = 1;
timelock1a = ft_selectdata(cfg, timelock1);
assert(isequal(size(timelock1a.avg), [1 5]));

cfg = [];
timelock2 = ft_appendtimelock(cfg, timelock1, timelock1, timelock1);

cfg = [];
cfg.channel = 1;
timelock2a = ft_selectdata(cfg, timelock2);
assert(isequal(size(timelock2a.trial), [3 1 5]));

cfg = [];
cfg.trials = [1 2];
timelock2b = ft_selectdata(cfg, timelock2);
assert(isequal(size(timelock2b.trial), [2 2 5]));

% The one that follows is a degenerate case. By selecting only one trial,
% the output is not really trial-based any more, but still contains one trial.
cfg = [];
cfg.trials = 1;
timelock2c = ft_selectdata(cfg, timelock2);
assert(isequal(size(timelock2c.trial), [1 2 5]));
% assert(isequal(size(timelock2c.trial), [2 5]));


%-------------------------------------
%generate data
data = [];
data.fsample = 1000;
data.cfg     = [];

nsmp  = 1000;
nchan = 80;
for k = 1:10
  data.trial{k} = randn(nchan,nsmp);
  data.time{k}  = ((1:nsmp)-1)./data.fsample;
end

% create grad-structure and add to data
grad.pnt  = randn(nchan,3);
grad.ori  = randn(nchan,3);
grad.tra  = eye(nchan);
for k = 1:nchan
  grad.label{k,1} = ['chan',num2str(k,'%03d')];
end
data.grad  = ft_datatype_sens(grad);
data.label = grad.label;
data.trialinfo = (1:10)';
data = ft_checkdata(data, 'hassampleinfo', 'yes');

%% this part of the script tests the functionality of ft_selectdata with respect
% to raw data.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this part of the script tests the functionality of ft_selectdata with respect
% to selecting the primary or a secondary dimord

freq = [];
freq.powspctrm = randn(3, 4, 5);
freq.dimord = 'chan_freq_time';
freq.crsspctrm = randn(3, 3, 4, 5);
freq.crsspctrmdimord = 'chan_chan_freq_time';
freq.label = {'1', '2', '3'};
freq.freq  = 1:4;
freq.time  = 1:5;

cfg = [];
freqpow = ft_selectdata(cfg, freq);

cfg = [];
cfg.parameter = 'powspctrm';
freqpow = ft_selectdata(cfg, freq);

cfg = [];
cfg.parameter = 'crsspctrm';
freqcrs = ft_selectdata(cfg, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this part of the script tests the functionality of ft_selectdata with respect
% to averaging over each dimension

% rpt_chan_freq_time

freq = [];
freq.dimord = 'rpt_chan_freq_time';
freq.label = {'1', '2', '3'};
freq.freq = 1:4;
freq.time = 1:5;
freq.powspctrm = randn(2, 3, 4, 5);

cfg = [];
cfg.avgoverrpt = 'yes';
cfg.keeprptdim = 'yes';
freq_avgoverrpt = ft_selectdata(cfg, freq);
cfg.keeprptdim = 'no';
freq_avgoverrpt = ft_selectdata(cfg, freq);

cfg = [];
cfg.avgoverchan = 'yes';
cfg.keepchandim = 'yes';
freq_avgoverchan = ft_selectdata(cfg, freq);
cfg.keepchandim = 'no';
freq_avgoverchan = ft_selectdata(cfg, freq);

cfg = [];
cfg.avgoverfreq = 'yes';
cfg.keepfreqdim = 'yes';
freq_avgoverfreq = ft_selectdata(cfg, freq);
cfg.keepfreqdim = 'no';
freq_avgoverfreq = ft_selectdata(cfg, freq);

cfg = [];
cfg.avgovertime = 'yes';
cfg.keeptimedim = 'yes';
freq_avgovertime = ft_selectdata(cfg, freq);
cfg.keeptimedim = 'no';
freq_avgovertime = ft_selectdata(cfg, freq);

cfg = [];
cfg.avgoverrpt  = 'yes';
cfg.avgoverchan = 'yes';
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'yes';
freq_avgoverall = ft_selectdata(cfg, freq);

% rpt_chan_time

timelock = [];
timelock.dimord = 'rpt_chan_time';
timelock.label = {'1', '2', '3'};
timelock.time = 1:4;
timelock.avg = randn(2, 3, 4);

cfg = [];
cfg.avgoverrpt = 'yes';
timelock_avgoverrpt = ft_selectdata(cfg, timelock);

cfg = [];
cfg.avgoverchan = 'yes';
timelock_avgoverchan = ft_selectdata(cfg, timelock);

cfg = [];
cfg.avgoverrpt = 'yes';
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
timelock_avgoverall = ft_selectdata(cfg, timelock);

% rpt_chan_time

cfg = [];
cfg.avgovertime = 'yes';
timelock_avgovertime = ft_selectdata(cfg, timelock);

timelock = [];
timelock.dimord = 'chan_time';
timelock.label = {'1', '2', '3'};
timelock.time = 1:4;
timelock.avg = randn(3, 4);

cfg = [];
cfg.avgoverchan = 'yes';
timelock_avgoverchan = ft_selectdata(cfg, timelock);

cfg = [];
cfg.avgovertime = 'yes';
timelock_avgovertime = ft_selectdata(cfg, timelock);

cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
timelock_avgoverall = ft_selectdata(cfg, timelock);

source = [];
source.dim = [10 11 12];
source.transform = eye(4);
source.avg.pow = rand(10*11*12,1);
source.inside = 1:660;
source.outside = 661:1320;

cfg = [];
cfg.avgoverpos = 'yes';
output = ft_selectdata(cfg, source);
assert(output.pos(1)==5.5);
assert(output.pos(2)==6.0);
assert(output.pos(3)==6.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this part of the script tests the functionality of ft_selectdata with respect
% to freqdata. it implements the (old) test_ft_selectdata_freqdata

% do spectral analysis
cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [2 100];
cfg.pad    = 1;
cfg.tapsmofrq = 3;
freq       = ft_freqanalysis(cfg, data);

cfg.output = 'pow';
cfg.keeptrials = 'yes';
freqp      = ft_freqanalysis(cfg, data);

cfg.output = 'powandcsd';
cfg.channelcmb = ft_channelcombination([data.label(1) {'all'};data.label(2) {'all'}], data.label);
freqc      = ft_freqanalysis(cfg, data);

cfg        = [];
cfg.method = 'mtmconvol';
cfg.foi    = [20:20:100]; % there are 10 repetitions, so let's use 5 frequencies
cfg.toi    = [0.4 0.5 0.6];
cfg.t_ftimwin = ones(1,numel(cfg.foi)).*0.2;
cfg.taper  = 'hanning';
cfg.output = 'pow';
cfg.keeptrials = 'yes';
freqtf     = ft_freqanalysis(cfg, data);


clear freq*

% make some dummy frequency structures
freq1.label = {'1' '2'};
freq1.freq  = 1:10;
freq1.time  = 1:5;
freq1.dimord = 'chan_freq_time';
freq1.powspctrm = randn(2,10,5);
freq1.cfg   = [];

cfg = [];
cfg.parameter = 'powspctrm';
freq2  = ft_appendfreq(cfg, freq1, freq1);
freq2  = rmfield(freq2, 'cfg');

clear freq*

freq3.label = {'1' '2'};
freq3.freq  = 1:10;
freq3.dimord = 'chan_freq';
freq3.powspctrm = randn(2,10);

cfg = [];
cfg.parameter = 'powspctrm';
freq4  = ft_appendfreq(cfg, freq3, freq3);
freq4  = rmfield(freq4, 'cfg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this part of the function tests the functionality of ft_selectdata with respect to timelock data

% create timelocked data
cfg = [];
cfg.keeptrials = 'yes';
tlck = ft_timelockanalysis(cfg, data);
cfg.covariance = 'yes';
tlckc = ft_timelockanalysis(cfg, data);
cfg.keeptrials = 'no';
tlckcavg = ft_timelockanalysis(cfg, data);
cfg.covariance = 'no';
tlckavg = ft_timelockanalysis(cfg, data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this part of the script tests the functionality of ft_selectdata with selections
% that are made into multiple fields present in the data

if false
  % this section does not yet work on 5 April 2014, so no point in testing
  freq = [];
  freq.powspctrm = randn(3, 4, 5);
  freq.dimord = 'chan_freq_time';
  freq.crsspctrm = randn(3, 3, 4, 5);
  freq.crsspctrmdimord = 'chan_chan_freq_time';
  freq.label = {'1', '2', '3'};
  freq.freq  = 1:4;
  freq.time  = 1:5;
  
  cfg = [];
  cfg.channel = {'1', '2'};
  cfg.foilim = [1 3];
  output = ft_selectdata(cfg, freq);
  assert(isfield(output, 'powspctrm'), 'field missing');
  assert(isfield(output, 'crsspctrm'), 'field missing');
  assert(size(output.powspctrm, 1)==2, 'incorrect size'); % chan
  assert(size(output.powspctrm, 2)==3, 'incorrect size'); % freq
  assert(size(output.crsspctrm, 1)==2, 'incorrect size'); % chan
  assert(size(output.crsspctrm, 2)==2, 'incorrect size'); % chan
  assert(size(output.crsspctrm, 3)==3, 'incorrect size'); % freq
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this part of the script tests the functionality of ft_selectdata with union/intersect on raw data

data = [];
data.trial{1} = randn(2,10);
data.time{1}  = -5:4;
data.trial{2} = randn(2,10);
data.time{2}  = 1:10;
data.label    = {'chan01';'chan02'};

cfg = [];
cfg.select = 'intersect';
datasel1 = ft_selectdata(cfg, data);

cfg.select = 'union';
ok = true;
try
  datasel2 = ft_selectdata(cfg, data);
catch
  ok = false;
end
assert(~ok, 'ft_selectdata with raw input and union selection should fail');
