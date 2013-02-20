function test_ft_selectdata

% TEST test_ft_selectdata
% TEST ft_selectdata ft_selectdata_old ft_selectdata_new ft_appendfreq

clear freq*

% make some dummy frequency structures
freq1.label = {'1' '2'};
freq1.freq  = 1:10;
freq1.time  = 1:5;
freq1.dimord = 'chan_freq_time';
freq1.powspctrm = randn(2,10,5);

cfg = [];
cfg.parameter = 'powspctrm';
freq2  = ft_appendfreq(cfg, freq1, freq1);
freq2  = rmfield(freq2, 'cfg');
freq2a = ft_selectdata(freq1, freq1, 'param', 'powspctrm'); % this should append the power spectrum
assert(isequal(freq2, freq2a));

freq4a = ft_selectdata(freq2, freq2, 'param', 'powspctrm');
assert(isequal(size(freq4a.powspctrm), [4 2 10 5]));

clear freq*

freq3.label = {'1' '2'};
freq3.freq  = 1:10;
freq3.dimord = 'chan_freq';
freq3.powspctrm = randn(2,10);

cfg = [];
cfg.parameter = 'powspctrm';
freq4  = ft_appendfreq(cfg, freq3, freq3);
freq4  = rmfield(freq4, 'cfg');
freq4a = ft_selectdata(freq3, freq3, 'param', 'powspctrm');  % this should append the power spectrum
assert(isequal(freq4, freq4a));

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


%% this part of the script tests the functionality of selectdata with respect
% to freqdata. it implements the (old) test_selectdata_freqdata 

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
  grad.label{k} = ['chan',num2str(k,'%03d')];
end
data.grad  = ft_datatype_sens(grad);
data.label = grad.label;

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
cfg.foi    = [10:10:100];
cfg.toi    = [0.4 0.5 0.6];
cfg.t_ftimwin = ones(1,numel(cfg.foi)).*0.2;
cfg.taper  = 'hanning';
cfg.output = 'pow';
cfg.keeptrials = 'yes';
freqtf     = ft_freqanalysis(cfg, data);

%select channels
fx1 = selectdata(freq,  'channel', data.label(5:10));
fp1 = selectdata(freqp, 'channel', data.label(5:10));
try
  fc1 = selectdata(freqc, 'channel', data.label(5:10)); % gives error
catch
  fprintf('selecting channels with csd in input does not work');
end
ftf1 = selectdata(freqtf, 'channel', data.label(5:10));

%select frequencies
fx2 = selectdata(freq,  'foilim', [10 40]);
fp2 = selectdata(freqp, 'foilim', [10 40]);
fc2 = selectdata(freqc, 'foilim', [10 40]);
ftf2 = selectdata(freqtf, 'foilim', [10 40]);

%select time
ftf2b = selectdata(freqtf, 'toilim', [0.5 0.6]);

%select trials
fx3 = selectdata(freq,  'rpt', 3:5);
fp3 = selectdata(freqp, 'rpt', 3:5);
fc3 = selectdata(freqc, 'rpt', 3:5);
ftf3 = selectdata(freqtf, 'rpt', 3:5);

%avgover channels
fx4 = selectdata(freq,  'avgoverchan', 'yes');
fp4 = selectdata(freqp, 'avgoverchan', 'yes');
fc4 = selectdata(freqc, 'avgoverchan', 'yes');
ftf4 = selectdata(freqtf, 'avgoverchan', 'yes');

%avgover frequencies
fx5 = selectdata(freq,  'avgoverfreq', 'yes');
fp5 = selectdata(freqp, 'avgoverfreq', 'yes');
fc5 = selectdata(freqc, 'avgoverfreq', 'yes');
ftf5 = selectdata(freqtf, 'avgoverchan', 'yes');

%avgover trials
fx6 = selectdata(freq,  'avgoverrpt', 'yes');
fp6 = selectdata(freqp, 'avgoverrpt', 'yes');
fc6 = selectdata(freqc, 'avgoverrpt', 'yes');
ftf6 = selectdata(freqtf, 'avgoverrpt', 'yes');

%leaveoneout
fx7 = selectdata(freq,  'jackknife', 'yes');
fp7 = selectdata(freqp, 'jackknife', 'yes');
fc7 = selectdata(freqc, 'jackknife', 'yes');
ftf7 = selectdata(freqtf, 'jackknife', 'yes');

