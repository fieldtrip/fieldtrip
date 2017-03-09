function test_bug103

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_singleplotER

freq.freq       = 1:1:100;
freq.powspctrm  = randn(size(freq.freq)).^2;
freq.label      = {'chan1'};
freq.dimord     = 'chan_freq';

cfg = [];
figure; ft_singleplotER(cfg, freq);

save /tmp/test_bug103.mat freq

cfg = [];
cfg.inputfile = '/tmp/test_bug103.mat';
figure; ft_singleplotER(cfg);
