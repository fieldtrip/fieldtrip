function test_bug893

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqstatistics statfun_indepsamplesZcoh

% function to confirm the bug and to test the fix.
% ft_freqstatistics fails when input data has 'rpttap' in the dimord

% create some data
freq               = [];
freq.label         = {'1';'2'};
freq.freq          = 1:10;
freq.fourierspctrm = randn(20,2,10)+randn(20,2,10).*1i;
freq.cumtapcnt     = 2*ones(10,1);
freq.dimord        = 'rpttap_chan_freq';

% create a configuration for the statistics
cfg           = [];
cfg.design    = [ones(1,10) ones(1,10)*2];
cfg.parameter = 'fourierspctrm';
cfg.statistic = 'ft_statfun_indepsamplesZcoh';
cfg.method    = 'montecarlo';
cfg.numrandomization = 1;
cfg.label     = freq.label;
stat          = ft_freqstatistics(cfg, freq);

% running this prior to the fix confirms the bug
