function test_bug1556

% MEM 1gb
% WALLTIME 00:03:04

% TEST test_bug1556
% TEST statfun_depsamplesF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1.powspctrm = [1     0     3     8     8    10     1     4     6     1]';
data1.dimord    = 'rpt_chan_freq_time';
data1.label     = {'cz'};
data1.freq      = 1;
data1.time      = 1;

data2.powspctrm = [1     8     6     1     6     8     7     4     9     0]';
data2.dimord    = 'rpt_chan_freq_time';
data2.label     = {'cz'};
data2.freq      = 1;
data2.time      = 1;

data3.powspctrm = [10     5     6     7    10     5    11     6    10     3]';
data3.dimord    = 'rpt_chan_freq_time';
data3.label     = {'cz'};
data3.freq      = 1;
data3.time      = 1;

cfg = [];
cfg.verbose  = 'off';
cfg.method   = 'analytic';
cfg.feedback = 'no';
cfg.alpha    = 5.0000e-02;
cfg.tail = 1;
cfg.statistic = 'depsamplesF';
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.design    = [ 1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2  3 3 3 3 3 3 3 3 3 3
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 ];

stats = ft_freqstatistics(cfg, data1, data2, data3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.verbose  = 'off';
cfg.method   = 'analytic';
cfg.feedback = 'no';
cfg.alpha    = 5.0000e-02;
cfg.tail = 1;
cfg.design    = [ 1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
                  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 ];

cfg.ivar      = 1;
cfg.statistic = 'indepsamplesT';
stats_indepsamplesT = ft_freqstatistics(cfg, data1, data2);
cfg.statistic = 'indepsamplesF';
stats_indepsamplesF = ft_freqstatistics(cfg, data1, data2);
cfg.uvar = 2;
cfg.statistic = 'depsamplesT';
stats_depsamplesT = ft_freqstatistics(cfg, data1, data2);
cfg.statistic = 'depsamplesF';
stats_depsamplesF = ft_freqstatistics(cfg, data1, data2);

stats_indepsamplesT.stat
stats_indepsamplesF.stat
stats_depsamplesT.stat
stats_depsamplesF.stat

dat = cat(1, data1.powspctrm, data2.powspctrm);
x = data1.powspctrm;
y = data2.powspctrm;

[p, anovatab, stats_f] = anova1([x, y], []);
[h, p, ci, stats_t2] = ttest2(x, y);
[h, p, ci, stats_t] = ttest(x-y);

if abs(stats_t2.tstat^2/anovatab{2,5} - 1) > 0.001
  error('t^2 is unequal to F');
end

