function test_issue2017

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_selectdata getdimord

freq = [];
freq.powspctrm = randn(10,1,1,1);
freq.powspctrm(:) = 1:10;
freq.dimord = 'rpt_chan_freq_time';
freq.label = {'1'};
freq.freq = 1;
freq.time = 1;

cfg = [];
cfg.trials = [1 3 5 7 9];
freqsel = ft_selectdata(cfg, freq);

assert(sum(freqsel.powspctrm(:))==sum(cfg.trials));
