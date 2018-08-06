function test_bug1729

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_multiplotTFR ft_singleplotTFR ft_plot_matrix

% reproduce the incorrect display
freq.label = {'chan'};
freq.time = [0 1 2];
freq.freq = [1:10 21:30];
freq.powspctrm = zeros(1, 20, 3);
freq.powspctrm(1,:,1) = freq.freq;
freq.powspctrm(1,:,2) = freq.freq;
freq.powspctrm(1,:,3) = freq.freq;
freq.dimord = 'chan_freq_time';

cfg = [];
ft_singleplotTFR(cfg, freq);

% try with another spacing, this also does not look
freq.label = {'chan'};
freq.time = [0 1 2];
freq.freq = logspace(log10(1),log10(30),20);
freq.powspctrm = zeros(1, 20, 3);
freq.powspctrm(1,:,1) = 1:20;
freq.powspctrm(1,:,2) = 1:20;
freq.powspctrm(1,:,3) = 1:20;
freq.dimord = 'chan_freq_time';

cfg = [];
ft_singleplotTFR(cfg, freq);

% the following causes a probably unrelated problem
% please also have a look at this
freq.label = {'chan'};
freq.time = 0;
freq.freq = [1:10 21:30];
freq.powspctrm = zeros(1, 20, 1);
freq.powspctrm(1,:,1) = freq.freq;
freq.dimord = 'chan_freq_time';

cfg = [];
ft_singleplotTFR(cfg, freq);
