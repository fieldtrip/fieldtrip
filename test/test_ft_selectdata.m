% TEST test_ft_selectdata ft_selectdata

clear freq*

% make some dummy frequency structures
freq1.label = {'1' '2'};
freq1.freq  = 1:10;
freq1.time  = 1:5;
freq1.dimord = 'chan_freq_time';
freq1.powspctrm = randn(2,10,5);

cfg = [];
cfg.parameter = 'powspctrm';
freq2 = ft_appendfreq(cfg, freq1, freq1);
freq2 = rmfield(freq2, 'cfg');
freq2a = ft_selectdata(freq2, freq2);
assert(identical(freq2, freq2a));
freq4a = ft_selectdata(freq2, freq2, 'param', 'powspctrm');
assert(identical(freq2, freq2a));

freq3.label = {'1' '2'};
freq3.freq  = 1:10;
freq3.dimord = 'chan_freq';
freq3.powspctrm = randn(2,10);

cfg = [];
cfg.parameter = 'powspctrm';
freq4 = ft_appendfreq(cfg, freq3, freq3);
freq4 = rmfield(freq4, 'cfg');
freq4a = ft_selectdata(freq3, freq3);
assert(identical(freq4, freq4a));
freq4a = ft_selectdata(freq3, freq3, 'param', 'powspctrm');
assert(identical(freq4, freq4a));

