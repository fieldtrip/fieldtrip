function test_ft_freqgrandaverage

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqgrandaverage

% this functions tests the new implementation of ft_freqgrandaverage. the
% new functionality includes the use of a cfg.parameter, and allows for
% 'chan' data to be averaged/combined within this function. it just tests
% whether it runs through for the most commonly used applications. it
% doesn't test the foilim/toilim selection

% create some data
freq1.label = {'chan1';'chan2'};
freq1.freq  = 1:10;
freq1.time  = 1:5;
freq1.dimord = 'chan_freq_time';
freq1.powspctrm = randn(2,10,5);
freq1.stat  = randn(2,10,5);

freq2 = freq1;
freq2.powspctrm = randn(2,10,5);
freq2.stat  = randn(2,10,5);

cfg = [];
grandavg1 = ft_freqgrandaverage(cfg, freq1, freq2);

cfg.parameter = 'stat';
grandavg2 = ft_freqgrandaverage(cfg, freq1, freq2);

cfg.parameter = {'powspctrm' 'stat'};
grandavg3 = ft_freqgrandaverage(cfg, freq1, freq2);

cfg.keepindividual = 'yes';
grandavg4 = ft_freqgrandaverage(cfg, freq1, freq2);

freq1 = rmfield(freq1, 'time');
freq2 = rmfield(freq2, 'time');
freq1.dimord = 'chan_freq';
freq2.dimord = 'chan_freq';
freq1.powspctrm = randn(2,10);
freq2.powspctrm = randn(2,10);
freq1.stat = randn(2,10);
freq2.stat = randn(2,10);


cfg.keepindividual = 'no';
grandavg5 = ft_freqgrandaverage(cfg, freq1, freq2);

cfg.keepindividual = 'yes';
grandavg6 = ft_freqgrandaverage(cfg, freq1, freq2);

freq1 = rmfield(freq1, 'freq');
freq2 = rmfield(freq2, 'freq');
freq1 = rmfield(freq1, 'powspctrm');
freq2 = rmfield(freq2, 'powspctrm');
freq1.dimord = 'chan';
freq2.dimord = 'chan';
freq1.stat = [1;2];
freq2.stat = [3;4];

cfg.parameter = 'stat';
cfg.keepindividual = 'no';
grandavg7 = ft_freqgrandaverage(cfg, freq1, freq2);

cfg.keepindividual = 'yes';
grandavg8 = ft_freqgrandaverage(cfg, freq1, freq2);
