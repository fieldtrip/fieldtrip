function test_bug1357

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_multiplotTFR ft_singleplotTFR ft_singleplotTFR ft_singleplotTFR

cfg = [];
cfg.layout = 'CTF151.lay';
lay = ft_prepare_layout(cfg);

label = ft_channelselection('MEG', lay.label);

nchan = length(label);
nfreq = 50;
ntime = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = [];
freq.label = label;
freq.freq = 1:nfreq;
freq.dimord = 'chan_freq';
freq.powspctrm = randn(nchan, nfreq);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.layout = lay;
ft_multiplotER(cfg, freq);
ft_singleplotER(cfg, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = [];
freq.label = label;
freq.freq = 1:nfreq;
freq.time = 1:ntime;
freq.dimord = 'chan_freq_time';
freq.powspctrm = randn(nchan, nfreq, ntime);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.layout = lay;
ft_multiplotTFR(cfg, freq);
ft_multiplotER(cfg, freq);
ft_singleplotTFR(cfg, freq);
ft_singleplotER(cfg, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = [];
freq.label = label;
freq.freq = 1:nfreq;
freq.time = 1:ntime;
freq.dimord = 'chan_time_freq'; % swap time and freq
freq.powspctrm = randn(nchan, ntime, nfreq);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.layout = lay;
ft_multiplotTFR(cfg, freq);
ft_multiplotER(cfg, freq);
ft_singleplotTFR(cfg, freq);
ft_singleplotER(cfg, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following one is what Giovanni explicitly reports

nchan = 1;
nfreq = 10;
ntime = 100;

freq = [];
freq.label = label(1:nchan);
freq.freq = 1:nfreq;
freq.time = 1:ntime;
freq.dimord = 'chan_freq_time';
freq.powspctrm = randn(nchan, nfreq, ntime);
ft_multiplotTFR(cfg, freq);
ft_multiplotER(cfg, freq);
ft_singleplotTFR(cfg, freq);
ft_singleplotER(cfg, freq);
