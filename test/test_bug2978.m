% function test_bug2978

% WALLTIME 00:20:00
% MEM 2GB

% TEST ft_multiplotER


dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare the data

cfg = [];
cfg.channel = 'MEG';
cfg.demean = 'yes';
cfg.dataset = dataset;
data = ft_preprocessing(cfg); % 266 trials

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
cfg.output = 'powandcsd'; % this implies keeptrials
powandcsd = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'coh';
coh1 = ft_connectivityanalysis(cfg, powandcsd);

% HACK, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3128
% coh1.dimord = 'chancmb_freq'; 

% I don't know how to get single trial (pseudo) estimates of coherence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, timelock

cfg = [];
cfg.layout = 'CTF151.lay';
figure; ft_multiplotER(cfg, timelock1);
figure; ft_multiplotER(cfg, data);
figure; ft_multiplotER(cfg, timelock2);


cfg = [];
cfg.layout = 'CTF151.lay';
cfg.trials = 255; % contains a jump
figure; ft_multiplotER(cfg, data);
figure; ft_multiplotER(cfg, timelock2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, power

cfg = [];
cfg.layout = 'CTF151.lay';
figure; ft_multiplotER(cfg, freq1);
figure; ft_multiplotER(cfg, freq2);

cfg = [];
cfg.layout = 'CTF151.lay';
cfg.trials = 255; % contains a jump
figure; ft_multiplotER(cfg, freq2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, timefreq

cfg = [];
cfg.layout = 'CTF151.lay';
figure; ft_multiplotER(cfg, timefreq1);
figure; ft_multiplotER(cfg, timefreq2);

cfg = [];
cfg.layout = 'CTF151.lay';
cfg.trials = 255; % contains a jump
figure; ft_multiplotER(cfg, timefreq2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, bivariate

cfg = [];
cfg.layout = 'CTF151.lay';
cfg.parameter = 'cohspctrm';
cfg.refchannel = 42;
figure; ft_multiplotER(cfg, coh1);



