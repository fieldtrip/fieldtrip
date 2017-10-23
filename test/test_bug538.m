function test_bug538

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqanalysis ft_connectivityanalysis ft_freqstatistics

% this script addresses bug 538, i.e. incompatibility between
% ft_freqstatistics, and connectivity data (with 'chan_chan' in dimord)

% create some data
for  k = 1:5
    trial1{k} = randn(2,100);
    trial2{k} = randn(2,100);
    trial3{k} = randn(2,100);
    trial4{k} = randn(2,100);
    time{k}   = (0:99)./100;
end
data1.trial = trial1;
data2.trial = trial2;
data3.trial = trial3;
data4.trial = trial4;
data1.time  = time;
data2.time  = time;
data3.time  = time;
data4.time  = time;
data1.label = {'1' '2'};
data2.label = {'1' '2'};
data3.label = {'1' '2'};
data4.label = {'1' '2'};

cfg = [];
cfg.method = 'mtmfft';
cfg.taper  = 'hanning';
cfg.output = 'fourier';
cfg.foilim = [0 40];
freq1 = ft_freqanalysis(cfg, data1);
freq2 = ft_freqanalysis(cfg, data2);
freq3 = ft_freqanalysis(cfg, data3);
freq4 = ft_freqanalysis(cfg, data4);

cfg = [];
cfg.method = 'plv';
plv1 = ft_connectivityanalysis(cfg, freq1);
plv2 = ft_connectivityanalysis(cfg, freq2);
plv3 = ft_connectivityanalysis(cfg, freq3);
plv4 = ft_connectivityanalysis(cfg, freq4);


cfg = [];
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.parameter = 'plvspctrm';
cfg.design    = [1 1 2 2;1 2 1 2];
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.method    = 'montecarlo';
cfg.numrandomization = 10;
stat = ft_freqstatistics(cfg, plv1, plv2, plv3, plv4);
