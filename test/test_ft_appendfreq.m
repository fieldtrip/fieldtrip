%TEST: test_ft_appendfreq ft_appendfreq

% make some dummy frequency structures
freq1.label = {'1';'2'};
freq1.freq  = 1:10;
freq1.time  = 1:5;
freq1.dimord = 'chan_freq_time';
freq1.powspctrm = randn(2,10,5);

cfg = [];

freq2         = freq1;
cfg.appenddim = 'rpt';
freqrpt       = ft_appendfreq(cfg, freq1, freq2);

freq2         = freq1;
freq2.label   = {'3';'4'};
cfg.appenddim = 'chan';
freqchan      = ft_appendfreq(cfg, freq1, freq2);

freq2         = freq1;
freq2.freq    = 11:20;
cfg.appenddim = 'freq';
freqfreq      = ft_appendfreq(cfg, freq1, freq2);

freq2         = freq1;
freq2.time    = 6:10;
cfg.appenddim = 'time';
freqtime      = ft_appendfreq(cfg, freq1, freq2);
