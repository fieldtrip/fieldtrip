function test_ft_selectdata

% MEM 1500mb
% WALLTIME 00:03:16

% TEST test_ft_selectdata
% TEST ft_selectdata ft_selectdata_old ft_selectdata_new ft_appendfreq

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

sel = [5 8 12 38];
cfg = [];cfg.channel = data.label(sel);
d1  = ft_selectdata(data, 'channel', data.label(sel));
d2  = ft_selectdata(cfg, data);
d2.cfg = [];
assert(isequal(d1,d2));

sel = [3 4 6 9];
cfg = [];cfg.trials = sel;
d3  = ft_selectdata(data, 'rpt', sel);
d4  = ft_selectdata(cfg, data);
d4.cfg = []; % ft_selectdata_old does not do anything with cfg
assert(isequal(d3,d4));

sel = [];
cfg = [];cfg.trials = sel;
d5  = ft_selectdata(data, 'rpt', sel);
d6  = ft_selectdata(cfg, data);
d6.cfg = []; % ft_selectdata_old does not do anything with cfg
assert(isequal(d5,d6));

sel = 'all';
cfg = [];cfg.trials = sel;
d7  = ft_selectdata(data, 'rpt', sel);
d8  = ft_selectdata(cfg, data);
d8.cfg = []; % ft_selectdata_old does not do anything with cfg
assert(isequal(d7,d8));
assert(isequal(d7,data));
assert(isequal(d8,data));

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
cfg.foi    = [10:10:100];
cfg.toi    = [0.4 0.5 0.6];
cfg.t_ftimwin = ones(1,numel(cfg.foi)).*0.2;
cfg.taper  = 'hanning';
cfg.output = 'pow';
cfg.keeptrials = 'yes';
freqtf     = ft_freqanalysis(cfg, data);

%% select channels, compare ft_selectdata_old with ft_selectdata_new and
% compare ft_selectdata_new with what would be expected

% make a selection of channels
sel  = data.label(5:10);
cfg  = [];cfg.channel = sel;
fx1  = ft_selectdata(freq,  'channel', data.label(5:10));
fx1b = ft_selectdata(cfg, freq);
fx1  = rmfield(fx1, 'cfg');
fx1b = rmfield(fx1b, 'cfg');
assert(isequal(fx1.fourierspctrm, freq.fourierspctrm(:,5:10,:)));
assert(isequal(fx1, fx1b));

fp1  = ft_selectdata(freqp, 'channel', data.label(5:10));
fp1b = ft_selectdata(cfg, freqp); 
fp1  = rmfield(fp1, 'cfg');
fp1b = rmfield(fp1b, 'cfg');
assert(isequal(fp1.powspctrm, freqp.powspctrm(:,5:10,:)));
assert(isequal(fp1, fp1b));

try
  fc1 = ft_selectdata(freqc, 'channel', data.label(5:10)); % gives error
catch
  fprintf('selecting channels with csd in input does not work');
end
ftf1  = ft_selectdata(freqtf, 'channel', data.label(5:10));
ftf1b = ft_selectdata(cfg, freqtf); 
ftf1  = rmfield(ftf1, 'cfg');
ftf1b = rmfield(ftf1b, 'cfg');
assert(isequal(ftf1.powspctrm, freqtf.powspctrm(:,5:10,:,:)));
assert(isequal(ftf1, ftf1b));

% select all channels
sel  = 'all';
cfg  = [];cfg.channel = sel;
fx1  = ft_selectdata(freq,  'channel', 'all');
fx1b = ft_selectdata(cfg, freq);
fx1  = rmfield(fx1, 'cfg');
fx1b = rmfield(fx1b, 'cfg');
assert(isequal(fx1.fourierspctrm, freq.fourierspctrm(:,1:80,:)));
assert(isequal(fx1, fx1b));

fp1  = ft_selectdata(freqp, 'channel', 'all');
fp1b = ft_selectdata(cfg, freqp); 
fp1  = rmfield(fp1, 'cfg');
fp1b = rmfield(fp1b, 'cfg');
assert(isequal(fp1.powspctrm, freqp.powspctrm(:,1:80,:)));
assert(isequal(fp1, fp1b));

try
  fc1 = ft_selectdata(freqc, 'channel', 'all'); % gives error
catch
  fprintf('selecting channels with csd in input does not work');
end
ftf1  = ft_selectdata(freqtf, 'channel', 'all');
ftf1b = ft_selectdata(cfg, freqtf); 
ftf1  = rmfield(ftf1, 'cfg');
ftf1b = rmfield(ftf1b, 'cfg');
assert(isequal(ftf1.powspctrm, freqtf.powspctrm(:,1:80,:,:)));
assert(isequal(ftf1, ftf1b));

% select no channels
sel  = [];
cfg  = [];cfg.channel = sel;
fx1  = ft_selectdata(freq,  'channel', []);
fx1b = ft_selectdata(cfg, freq);
fx1  = rmfield(fx1, 'cfg');
fx1b = rmfield(fx1b, 'cfg');
%assert(isequal(fx1.fourierspctrm, freq.fourierspctrm(:,5:10,:)));
assert(isequal(fx1, fx1b));

fp1  = ft_selectdata(freqp, 'channel', []);
fp1b = ft_selectdata(cfg, freqp); 
fp1  = rmfield(fp1, 'cfg');
fp1b = rmfield(fp1b, 'cfg');
%assert(isequal(fp1.powspctrm, freqp.powspctrm(:,5:10,:)));
assert(isequal(fp1, fp1b));

try
  fc1 = ft_selectdata(freqc, 'channel', []); % gives error
catch
  fprintf('selecting channels with csd in input does not work');
end
ftf1  = ft_selectdata(freqtf, 'channel', []);
ftf1b = ft_selectdata(cfg, freqtf); 
ftf1  = rmfield(ftf1, 'cfg');
ftf1b = rmfield(ftf1b, 'cfg');
%assert(isequal(ftf1.powspctrm, freqtf.powspctrm(:,5:10,:,:)));
assert(isequal(ftf1, ftf1b));

%% select frequencies

% subselection
sel  = [10 40];
cfg  = []; cfg.frequency = sel;
fx2  = ft_selectdata(freq,  'foilim', sel);
fx2b = ft_selectdata(cfg, freq);
fx2  = rmfield(fx2, 'cfg');
fx2b = rmfield(fx2b, 'cfg');
assert(isequal(fx2.fourierspctrm, freq.fourierspctrm(:,:,9:39)));
assert(isequal(fx2, fx2b));

fp2  = ft_selectdata(freqp, 'foilim', sel);
fp2b = ft_selectdata(cfg, freqp); 
fp2  = rmfield(fp2, 'cfg');
fp2b = rmfield(fp2b, 'cfg');
assert(isequal(fp2.powspctrm, freqp.powspctrm(:,:,9:39)));
assert(isequal(fp2, fp2b));

fc2 = ft_selectdata(freqc, 'foilim', sel);

ftf2  = ft_selectdata(freqtf, 'foilim', sel);
ftf2b = ft_selectdata(cfg, freqtf); 
ftf2  = rmfield(ftf2, 'cfg');
ftf2b = rmfield(ftf2b, 'cfg');
assert(isequal(ftf2.powspctrm, freqtf.powspctrm(:,:,1:4,:)));
assert(isequal(ftf2, ftf2b));

% all frequencies
sel  = 'all';
cfg  = []; cfg.frequency = sel;
fx2  = ft_selectdata(freq,  'foilim', sel);
fx2b = ft_selectdata(cfg, freq);
fx2  = rmfield(fx2, 'cfg');
fx2b = rmfield(fx2b, 'cfg');
assert(isequal(fx2.fourierspctrm, freq.fourierspctrm));
assert(isequal(fx2, fx2b));

fp2  = ft_selectdata(freqp, 'foilim', sel);
fp2b = ft_selectdata(cfg, freqp); 
fp2  = rmfield(fp2, 'cfg');
fp2b = rmfield(fp2b, 'cfg');
assert(isequal(fp2.powspctrm, freqp.powspctrm));
assert(isequal(fp2, fp2b));

fc2 = ft_selectdata(freqc, 'foilim', sel);

ftf2  = ft_selectdata(freqtf, 'foilim', sel);
ftf2b = ft_selectdata(cfg, freqtf); 
ftf2  = rmfield(ftf2, 'cfg');
ftf2b = rmfield(ftf2b, 'cfg');
assert(isequal(ftf2.powspctrm, freqtf.powspctrm));
assert(isequal(ftf2, ftf2b));

% no frequencies
sel  = [];
cfg  = []; cfg.frequency = sel;
fx2  = ft_selectdata(freq,  'foilim', sel);
fx2b = ft_selectdata(cfg, freq);
fx2  = rmfield(fx2, 'cfg');
fx2b = rmfield(fx2b, 'cfg');
%assert(isequal(fx2.fourierspctrm, freq.fourierspctrm(:,:,9:39)));
assert(isequal(fx2, fx2b));

fp2  = ft_selectdata(freqp, 'foilim', sel);
fp2b = ft_selectdata(cfg, freqp); 
fp2  = rmfield(fp2, 'cfg');
fp2b = rmfield(fp2b, 'cfg');
%assert(isequal(fp2.powspctrm, freqp.powspctrm(:,:,9:39)));
assert(isequal(fp2, fp2b));

fc2 = ft_selectdata(freqc, 'foilim', sel);

ftf2  = ft_selectdata(freqtf, 'foilim', sel);
ftf2b = ft_selectdata(cfg, freqtf); 
ftf2  = rmfield(ftf2, 'cfg');
ftf2b = rmfield(ftf2b, 'cfg');
%assert(isequal(ftf2.powspctrm, freqtf.powspctrm(:,:,1:4,:)));
assert(isequal(ftf2, ftf2b));

%% select time

% subselection
sel = [0.5 0.6]; selindx = [nearest(freqtf.time,0.5) nearest(freqtf.time,0.6)];
cfg = []; cfg.latency = sel;
ftf2  = ft_selectdata(freqtf, 'toilim', sel);
ftf2b = ft_selectdata(cfg, freqtf);
ftf2  = rmfield(ftf2, 'cfg');
ftf2b = rmfield(ftf2b, 'cfg');
assert(isequal(ftf2.powspctrm, freqtf.powspctrm(:,:,:,selindx(1):selindx(2))));
assert(isequal(ftf2, ftf2b));

% all
sel = 'all';
cfg = []; cfg.latency = sel;
ftf2  = ft_selectdata(freqtf, 'toilim', sel);
ftf2b = ft_selectdata(cfg, freqtf);
ftf2  = rmfield(ftf2, 'cfg');
ftf2b = rmfield(ftf2b, 'cfg');
assert(isequal(ftf2.powspctrm, freqtf.powspctrm));
assert(isequal(ftf2, ftf2b));

% nothing
sel = [];
cfg = []; cfg.latency = sel;
ftf2  = ft_selectdata(freqtf, 'toilim', sel);
ftf2b = ft_selectdata(cfg, freqtf);
ftf2  = rmfield(ftf2, 'cfg');
ftf2b = rmfield(ftf2b, 'cfg');
%assert(isequal(ftf2.powspctrm, freqtf.powspctrm(:,:,:,selindx(1):selindx(2))));
assert(isequal(ftf2, ftf2b));

%% select trials

% do a subselection
sel  = 3:5;
cfg  = []; cfg.trials = sel;
fx3  = ft_selectdata(freq,  'rpt', sel);
fp3  = ft_selectdata(freqp, 'rpt', sel);
fc3  = ft_selectdata(freqc, 'rpt', sel);
ftf3 = ft_selectdata(freqtf, 'rpt', sel);
fx3n  = ft_selectdata(cfg, freq);
fp3n  = ft_selectdata(cfg, freqp);
fc3n  = ft_selectdata(cfg, freqc);
ftf3n = ft_selectdata(cfg, freqtf);

fx3  = rmfield(fx3, 'cfg');
fx3n = rmfield(fx3n, 'cfg');
assert(isequal(fx3, fx3n));
fp3  = rmfield(fp3, 'cfg');
fp3n = rmfield(fp3n, 'cfg');
assert(isequal(fp3, fp3n));
fc3  = rmfield(fc3, 'cfg');
fc3n = rmfield(fc3n, 'cfg');
try, 
  assert(isequal(fc3, fc3n)); 
catch
  warning('assertion failed, because ft_selectdata_new cannot deal with crsspctrm in input yet');
end
ftf3  = rmfield(ftf3, 'cfg');
ftf3n = rmfield(ftf3n, 'cfg');
assert(isequal(ftf3, ftf3n));

% do an empty selection
sel  = [];
cfg  = []; cfg.trials = sel;
fx3  = ft_selectdata(freq,  'rpt', sel);
fp3  = ft_selectdata(freqp, 'rpt', sel);
fc3  = ft_selectdata(freqc, 'rpt', sel);
ftf3 = ft_selectdata(freqtf, 'rpt', sel);
fx3n  = ft_selectdata(cfg, freq);
fp3n  = ft_selectdata(cfg, freqp);
fc3n  = ft_selectdata(cfg, freqc);
ftf3n = ft_selectdata(cfg, freqtf);

fx3  = rmfield(fx3, 'cfg');
fx3n = rmfield(fx3n, 'cfg');
assert(isequal(fx3, fx3n));
fp3  = rmfield(fp3, 'cfg');
fp3n = rmfield(fp3n, 'cfg');
assert(isequal(fp3, fp3n));
fc3  = rmfield(fc3, 'cfg');
fc3n = rmfield(fc3n, 'cfg');
try, 
  assert(isequal(fc3, fc3n)); 
catch
  warning('assertion failed, because ft_selectdata_new cannot deal with crsspctrm in input yet');
end
ftf3  = rmfield(ftf3, 'cfg');
ftf3n = rmfield(ftf3n, 'cfg');
assert(isequal(ftf3, ftf3n));

% select all
sel  = 'all';
cfg  = []; cfg.trials = sel;
fx3  = ft_selectdata(freq,  'rpt', sel);
fp3  = ft_selectdata(freqp, 'rpt', sel);
fc3  = ft_selectdata(freqc, 'rpt', sel);
ftf3 = ft_selectdata(freqtf, 'rpt', sel);
fx3n  = ft_selectdata(cfg, freq);
fp3n  = ft_selectdata(cfg, freqp);
fc3n  = ft_selectdata(cfg, freqc);
ftf3n = ft_selectdata(cfg, freqtf);

fx3  = rmfield(fx3, 'cfg');
fx3n = rmfield(fx3n, 'cfg');
assert(isequal(fx3, fx3n));
fp3  = rmfield(fp3, 'cfg');
fp3n = rmfield(fp3n, 'cfg');
assert(isequal(fp3, fp3n));
fc3  = rmfield(fc3, 'cfg');
fc3n = rmfield(fc3n, 'cfg');
try, 
  assert(isequal(fc3, fc3n)); 
catch
  warning('assertion failed, because ft_selectdata_new cannot deal with crsspctrm in input yet');
end
ftf3  = rmfield(ftf3, 'cfg');
ftf3n = rmfield(ftf3n, 'cfg');
assert(isequal(ftf3, ftf3n));

%% avgover channels
fx4  = ft_selectdata(freq,   'avgoverchan', 'yes');
fp4  = ft_selectdata(freqp,  'avgoverchan', 'yes');
fc4  = ft_selectdata(freqc,  'avgoverchan', 'yes');
ftf4 = ft_selectdata(freqtf, 'avgoverchan', 'yes');

% assessing label after averaging: see bug 2191 -> this seems OK
cfg             = [];
cfg.avgoverchan = 'yes';
fx42  = ft_selectdata(cfg,freq);
fp42  = ft_selectdata(cfg,freqp);
fc42  = ft_selectdata(cfg,freqc);
ftf42 = ft_selectdata(cfg,freqtf);

if ~strcmp(fx4.label{:},fx42.label{:});error('mismatch on label field');end
if ~strcmp(fp4.label{:},fp42.label{:});error('mismatch on label field');end
if ~strcmp(fc4.label{:},fc42.label{:});error('mismatch on label field');end
if ~strcmp(ftf4.label{:},ftf42.label{:});error('mismatch on label field');end

fx4  = rmfield(fx4,  'cfg');
fx42 = rmfield(fx42, 'cfg');
assert(isequal(fx4, fx42));
fp4  = rmfield(fp4,  'cfg');
fp42 = rmfield(fp42, 'cfg');
assert(isequal(fp4, fp42));
fc4  = rmfield(fc4,  'cfg');
fc42 = rmfield(fc42, 'cfg');
try, 
  assert(isequal(fc4, fc42));
catch
  warning('ft_selectdata_new cannot work yet with crsspctrm');
end
ftf4  = rmfield(ftf4,  'cfg');
ftf42 = rmfield(ftf42, 'cfg');
assert(isequal(ftf4, ftf42));

%% avgover frequencies
fx5 = ft_selectdata(freq,  'avgoverfreq', 'yes');
fp5 = ft_selectdata(freqp, 'avgoverfreq', 'yes');
fc5 = ft_selectdata(freqc, 'avgoverfreq', 'yes');
ftf5 = ft_selectdata(freqtf, 'avgoverfreq', 'yes');

cfg             = [];
cfg.avgoverfreq = 'yes';
fx52  = ft_selectdata(cfg,freq);
fp52  = ft_selectdata(cfg,freqp);
fc52  = ft_selectdata(cfg,freqc);
ftf52 = ft_selectdata(cfg,freqtf);

fx5  = rmfield(fx5,  'cfg');
fx52 = rmfield(fx52, 'cfg');
fx5.dimord = fx52.dimord; % ft_selectdata_new displays the correct dimord, don't spend time on fixing this for ft_selectdata_old 
assert(isequal(fx5, fx52));
fp5  = rmfield(fp5,  'cfg');
fp52 = rmfield(fp52, 'cfg');
fp5.dimord = fp52.dimord; % ft_selectdata_new displays the correct dimord, don't spend time on fixing this for ft_selectdata_old 
assert(isequal(fp5, fp52));
fc5  = rmfield(fc5,  'cfg');
fc52 = rmfield(fc52, 'cfg');
fc5.dimord = fc52.dimord; % ft_selectdata_new displays the correct dimord, don't spend time on fixing this for ft_selectdata_old 
try, 
  assert(isequal(fc5, fc52));
catch
  warning('ft_selectdata_new cannot work yet with crsspctrm');
end
ftf5  = rmfield(ftf5,  'cfg');
ftf52 = rmfield(ftf52, 'cfg');
ftf5.dimord = ftf52.dimord; % ft_selectdata_new displays the correct dimord, don't spend time on fixing this for ft_selectdata_old 
assert(isequal(ftf5, ftf52));

%% avgover trials
fx6 = ft_selectdata(freq,  'avgoverrpt', 'yes');
fp6 = ft_selectdata(freqp, 'avgoverrpt', 'yes');
fc6 = ft_selectdata(freqc, 'avgoverrpt', 'yes');
ftf6 = ft_selectdata(freqtf, 'avgoverrpt', 'yes');

cfg            = [];
cfg.avgoverrpt = 'yes';
fx62  = ft_selectdata(cfg,freq);
fp62  = ft_selectdata(cfg,freqp);
fc62  = ft_selectdata(cfg,freqc);
ftf62 = ft_selectdata(cfg,freqtf);

fx6  = rmfield(fx6,  'cfg');
fx62 = rmfield(fx62, 'cfg');
%fx6.dimord = fx62.dimord; % ft_selectdata_new displays the correct dimord, don't spend time on fixing this for ft_selectdata_old 
fx6 = rmfield(fx6, {'cumtapcnt', 'cumsumcnt', 'trialinfo'}); % ft_selectdata_old does something with these fields, but ft_selectdata_new removes them. I would say that appropriate behavior is remove them
assert(isequal(fx6, fx62));
assert(isequal(fx6.fourierspctrm, squeeze(mean(freq.fourierspctrm))));
fp6  = rmfield(fp6,  'cfg');
fp62 = rmfield(fp62, 'cfg');
%fp6.dimord = fp62.dimord; % ft_selectdata_new displays the correct dimord, don't spend time on fixing this for ft_selectdata_old 
fp6 = rmfield(fp6, {'cumtapcnt', 'cumsumcnt', 'trialinfo'}); % ft_selectdata_old does something with these fields, but ft_selectdata_new removes them. I would say that appropriate behavior is remove them
assert(isequal(fp6, fp62));
assert(isequal(fp6.powspctrm, squeeze(mean(freqp.powspctrm))));
fc6  = rmfield(fc6,  'cfg');
fc62 = rmfield(fc62, 'cfg');
%fc6.dimord = fc62.dimord; % ft_selectdata_new displays the correct dimord, don't spend time on fixing this for ft_selectdata_old 
fc6 = rmfield(fc6, {'cumtapcnt', 'cumsumcnt', 'trialinfo'}); % ft_selectdata_old does something with these fields, but ft_selectdata_new removes them. I would say that appropriate behavior is remove them
try, 
  assert(isequal(fc6, fc62));
catch
  warning('ft_selectdata_new cannot work yet with crsspctrm');
end
ftf6  = rmfield(ftf6,  'cfg');
ftf62 = rmfield(ftf62, 'cfg');
ftf6.dimord = ftf62.dimord; % ft_selectdata_new displays the correct dimord, don't spend time on fixing this for ft_selectdata_old 
ftf6 = rmfield(ftf6, {'cumtapcnt', 'trialinfo'}); % ft_selectdata_old does something with these fields, but ft_selectdata_new removes them. I would say that appropriate behavior is remove them
assert(isequal(ftf6, ftf62));
assert(isequal(ftf6.powspctrm, squeeze(mean(freqtf.powspctrm))));

%% leaveoneout
% fx7 = ft_selectdata(freq,  'jackknife', 'yes'); %FAILS due to 'rpttap'
fp7 = ft_selectdata(freqp, 'jackknife', 'yes');
fc7 = ft_selectdata(freqc, 'jackknife', 'yes');
ftf7 = ft_selectdata(freqtf, 'jackknife', 'yes');

%% this part tests the functionality of ft_appendfreq

whos
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


