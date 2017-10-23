function test_bug1514

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_spike_select ft_selectdata

% the following structure corresponds to the one explained in
% ft_datatype_spike

spike = [];
spike.label = {'unit1' 'unit2' 'unit3'};
spike.timestamp = {randn(1,504) randn(1,50) randn(1,101)};
spike.time = {randn(1,504) randn(1,50) randn(1,101)};
spike.trial = {ceil(rand(1,504)*100) ceil(rand(1,50)*100) ceil(randn(1,101)*100)};
spike.trialtime = 3*randn(100,2);
spike.waveform = {randn(1,32,504) randn(1,32,50) randn(1,32,101)};
spike.waveformdimord = '{chan}_lead_time_spike';
spike.fourierspctrm = {randn(504,2,20), randn(50,2,20), randn(101,2,20)};
spike.fourierspctrmdimord = '{chan}_spike_lfplabel_freq';
spike.lfplabel = {'lfpchan1', 'lfpchan2'};
spike.freq = 1:20;


cfg = [];
cfg.trials = 1:50;
output = ft_selectdata(cfg, spike);

cfg = [];
cfg.channel = 2;
output = ft_selectdata(cfg, spike);

cfg = [];
cfg.latency = 'prestim';
output = ft_selectdata(cfg, spike);

cfg = [];
cfg.trials = 1:50;
cfg.channel = 2;
cfg.latency = 'prestim';
output = ft_selectdata(cfg, spike);

assert(all(output.time{1}<=0));
assert(all(output.trial{1}>=1 & output.trial{1}<=50));
assert(isequal(output.label, {'unit2'}));



