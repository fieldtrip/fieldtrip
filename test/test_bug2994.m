function test_bug2994

% TEST ft_freqstatistics
% TEST ft_freqgrandaverage
% TEST ft_selectdata

% WALLTIME 00:10:00
% MEM 250mb

% Bug 2994 has been reported by Niels. Apparently once upon a time it was
% possible to call ft_freqstatistics with ft_freqgrandaveraged data in the
% input, where the dimord is 'subj_freq_time' (NOTE the lacking chan
% dimord) with the data containing apparently an average across channels.
% The input data contains a label field. ft_freqstatistics crashes when the
% low level function tries to cluster across space.

% Let's start this form the beginning, create a multisubject
% freq-structure, with data averaged across channels, and see whether it
% pulls through

tmp.freq   = [1 2 3];
tmp.time   = [1 2 3 4 5];
tmp.label  = {'a';'b'};
tmp.dimord = 'chan_freq_time';


freq = cell(1,10);
freq0 = cell(1,10);
for k = 1:10
  tmp.powspctrm = rand(2,3,5);
  freq{k} = tmp;
  tmp.powspctrm(:) = 0;
  freq0{k} = tmp;  
end

%% ------------------------------------------------------------------------
% here, the following is done:
% -ft_freqgrandaverage (keepindividual = 'yes')
% -ft_selectdata (avgoverchan = 'yes')
% -ft_freqstatistics

cfg1 = [];
cfg1.keepindividual = 'yes';

cfg2 = [];
cfg2.avgoverchan = 'yes';

cfg3 = [];
cfg3.method    = 'montecarlo';
cfg3.statistic = 'ft_statfun_depsamplesT';
cfg3.design    = [ones(1,10) ones(1,10)*2;1:10 1:10];
cfg3.numrandomization = 100;
cfg3.correctm  = 'cluster';
cfg3.ivar      = 1;
cfg3.uvar      = 2;

freqall  = ft_selectdata(cfg2, ft_freqgrandaverage(cfg1, freq{:}));
freqall0 = ft_selectdata(cfg2, ft_freqgrandaverage(cfg1, freq0{:}));
assert(isequal(size(freqall.powspctrm),[10 1 3 5])); % there should be a singleton channel dimension
stat1    = ft_freqstatistics(cfg3, freqall, freqall0);

%% ------------------------------------------------------------------------
% here, the following is done:
% -ft_freqgrandaverage (keepindividual = 'yes')
% -ft_freqstatistics (avgoverchan = 'yes')

cfg1 = [];
cfg1.keepindividual = 'yes';

cfg3 = [];
cfg3.method    = 'montecarlo';
cfg3.statistic = 'ft_statfun_depsamplesT';
cfg3.avgoverchan = 'yes';
cfg3.design    = [ones(1,10) ones(1,10)*2;1:10 1:10];
cfg3.numrandomization = 100;
cfg3.correctm  = 'cluster';
cfg3.ivar      = 1;
cfg3.uvar      = 2;

freqall  = ft_freqgrandaverage(cfg1, freq{:});
freqall0 = ft_freqgrandaverage(cfg1, freq0{:});
assert(isequal(size(freqall.powspctrm),[10 2 3 5])); % there should be two channels
stat2    = ft_freqstatistics(cfg3, freqall, freqall0);

%% ------------------------------------------------------------------------
% here, the following is done:
% - NO ft_freqgrandaverage
% -ft_freqstatistics (avgoverchan = 'yes')

cfg3 = [];
cfg3.method    = 'montecarlo';
cfg3.statistic = 'ft_statfun_depsamplesT';
cfg3.avgoverchan = 'yes';
cfg3.design    = [ones(1,10) ones(1,10)*2;1:10 1:10];
cfg3.numrandomization = 100;
cfg3.correctm  = 'cluster';
cfg3.ivar      = 1;
cfg3.uvar      = 2;
cfg3.avgoverchan = 'yes';

stat3    = ft_freqstatistics(cfg3, freq{:}, freq0{:});

%% ------------------------------------------------------------------------
% here, the following is done:
% -ft_appendfreq
% -ft_freqstatistics (avgoverchan = 'yes')

cfg1 = [];
cfg1.parameter = 'powspctrm';
cfg1.appenddim = 'rpt';

cfg3 = [];
cfg3.method    = 'montecarlo';
cfg3.statistic = 'ft_statfun_depsamplesT';
cfg3.avgoverchan = 'yes';
cfg3.design    = [ones(1,10) ones(1,10)*2;1:10 1:10];
cfg3.numrandomization = 100;
cfg3.correctm  = 'cluster';
cfg3.ivar      = 1;
cfg3.uvar      = 2;
cfg3.avgoverchan = 'yes';

freqall  = ft_appendfreq(cfg1, freq{:});
freqall0 = ft_appendfreq(cfg1, freq0{:});

stat4    = ft_freqstatistics(cfg3, freqall, freqall0);

assert(isequal(stat1.stat,stat2.stat));
assert(isequal(stat1.stat,stat3.stat));
assert(isequal(stat1.stat,stat4.stat));
assert(isequal(stat2.stat,stat3.stat));
assert(isequal(stat2.stat,stat4.stat));
assert(isequal(stat3.stat,stat4.stat));

% the conclusion so far is that with the current state of the code, all
% seems to work fine.



