function test_bug3345

% MEM 5GB
% WALLTIME 00:20:00
% DEPENDENCY ft_multiplotER ft_singleplotER ft_topoplotER ft_multiplotTFR ft_singleplotTFR ft_topoplotTFR topolot_common bivariate_common chanscale_common

% See also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3296 
% which is on cfg.trials in the plotting functions

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare the data

cfg = [];
cfg.channel = 'MEG';
cfg.demean = 'yes';
cfg.dataset = dataset;
data = ft_preprocessing(cfg); % 266 trials

cfg = [];
cfg.trials = 1:10;
data = ft_selectdata(cfg, data); % only 10 trials to speed it all up
data.trial{10}(:,450:end) = data.trial{10}(:,450:end) + 0.2 * 1e-11; % insert a jump in trial 10

cfg = [];
timelock1 = ft_timelockanalysis(cfg, data);
cfg.keeptrials = 'yes';
timelock2 = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
freq1 = ft_freqanalysis(cfg, data);
cfg.keeptrials = 'yes';
freq2 = ft_freqanalysis(cfg, data);

cfg = [];
cfg.foilim = [1 20];
cfg.method = 'wavelet';
cfg.toi = -1:0.100:2;
timefreq1 = ft_freqanalysis(cfg, data);
cfg.keeptrials = 'yes';
timefreq2 = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'fourier'; % this implies keeptrials
cfg.foi = 1:10;
fourier = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'coh';
coh1 = ft_connectivityanalysis(cfg, fourier);

% HACK
coh2 = coh1;
coh2.cohspctrm = repmat(reshape(coh2.cohspctrm, 1, 151, 151, 10), 10, 1, 1, 1);
coh2.dimord = 'rpt_chan_chan_freq';

cfg = [];
cfg.method = 'wavelet';
cfg.output = 'fourier'; % this implies keeptrials
cfg.foi = 1:10;
cfg.toi = 0:0.1:2;
fourier = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'coh';
cohft1 = ft_connectivityanalysis(cfg, fourier);

% HACK
cohft2 = cohft1;
cohft2.cohspctrm = repmat(reshape(cohft2.cohspctrm, 1, 151, 151, 10, 21), 10, 1, 1, 1);
cohft2.dimord = 'rpt_chan_chan_freq_time';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, timelock
close all

plotoption = {@ft_multiplotER, @ft_singleplotER, @ft_topoplotER};
for xxx=1:3
  ft_xxxplot = plotoption{xxx};
  
  cfg = [];
  cfg.layout = 'CTF151.lay';
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, timelock1);
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, timelock1); % same as previous
  try
    errordetected = true;
    cfg.trials = 10;
    figure; ft_xxxplot(cfg, timelock1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, timelock1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, timelock2); % no jump
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, timelock2); % small jump (average)
  cfg.trials = 10;
  figure; ft_xxxplot(cfg, timelock2); % big jump
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, timelock2); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, data); % no jump
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, data); % small jump (average)
  cfg.trials = 10;
  figure; ft_xxxplot(cfg, data); % big jump
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, data); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
end % plotoption

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, power
close all

plotoption = {@ft_multiplotER, @ft_singleplotER, @ft_topoplotER};
for xxx=1:3
  ft_xxxplot = plotoption{xxx};
  
  cfg = [];
  cfg.layout = 'CTF151.lay';
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, freq1);
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, freq1); % same as previous
  try
    errordetected = true;
    cfg.trials = 10;
    figure; ft_xxxplot(cfg, freq1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, freq1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, freq2); % no jump
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, freq2); % small jump (average)
  cfg.trials = 10;
  figure; ft_xxxplot(cfg, freq2); % big jump
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, freq2); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
end % plotoption

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, timefreq
close all

plotoption = {@ft_multiplotTFR, @ft_singleplotTFR, @ft_topoplotTFR};
for xxx=1:3
  ft_xxxplot = plotoption{xxx};
  
  cfg = [];
  cfg.layout = 'CTF151.lay';
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, timefreq1);
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, timefreq1); % same as previous
  try
    errordetected = true;
    cfg.trials = 10;
    figure; ft_xxxplot(cfg, timefreq1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, timefreq1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, timefreq2); % no jump
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, timefreq2); % small jump (average)
  cfg.trials = 10;
  figure; ft_xxxplot(cfg, timefreq2); % big jump
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, timefreq2); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
end % plotoption

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, timefreq, plot the timecourse for a single frequency

close all

plotoption = {@ft_multiplotER, @ft_singleplotER};
for xxx=1:2
  ft_xxxplot = plotoption{xxx};
  
  cfg = [];
  cfg.layout = 'CTF151.lay';
  cfg.frequency = 10;
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, timefreq1);
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, timefreq1); % same as previous
  try
    errordetected = true;
    cfg.trials = 10;
    figure; ft_xxxplot(cfg, timefreq1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, timefreq1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, timefreq2); % no jump
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, timefreq2); % small jump (average)
  cfg.trials = 10;
  figure; ft_xxxplot(cfg, timefreq2); % big jump
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, timefreq2); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
end % plotoption

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, bivariate freq
close all

plotoption = {@ft_multiplotER, @ft_singleplotER, @ft_topoplotER};
for xxx=1:3
  ft_xxxplot = plotoption{xxx};
  
  cfg = [];
  cfg.layout = 'CTF151.lay';
  cfg.refchannel = 'gui';
  cfg.parameter = 'cohspctrm';
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, coh1);
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, coh1); % same as previous
  try
    errordetected = true;
    cfg.trials = 10;
    figure; ft_xxxplot(cfg, coh1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, coh1); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
  cfg.trials = 1;
  figure; ft_xxxplot(cfg, coh2); % no jump
  cfg.trials = 'all';
  figure; ft_xxxplot(cfg, coh2); % small jump (average)
  cfg.trials = 10;
  figure; ft_xxxplot(cfg, coh2); % big jump
  try
    errordetected = true;
    cfg.trials = [];
    figure; ft_xxxplot(cfg, coh2); % error
    errordetected = false;
  end
  assert(errordetected, 'this should have given an error');
  
end % plotoption

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, bivariate timefreq

plotoption = {@ft_multiplotTFR, @ft_singleplotTFR, @ft_topoplotTFR};
for xxx=1:3
  ft_xxxplot = plotoption{xxx};
  
  % FIXME
  
end % plotoption

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try time-freq and freq-time for the axes

plotoption = {@ft_multiplotTFR, @ft_singleplotTFR, @ft_topoplotTFR};
for xxx=1:3
  ft_xxxplot = plotoption{xxx};
  
  cfg = [];
  cfg.layout = 'CTF151.lay';
  cfg.parameter = 'powspctrm';
  cfg.xparam = 'time';
  cfg.yparam = 'freq';
  figure; ft_xxxplot(cfg, timefreq1);
  
  cfg.xparam = 'freq';
  cfg.yparam = 'time';
  figure; ft_xxxplot(cfg, timefreq1);
  
end
