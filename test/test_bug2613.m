function test_bug2613

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_freqanalysis ft_checkdata
% DATA private

load(dccnpath('/project/3031000.02/test/bug2613.mat'));

% cfg.trials is incorrect in the original as the selection runs all the way up to trial 240 (there are in fact 182 trials)
% this fixes it and also speeds up the test
cfg.trials = 1:10;

freq = ft_freqanalysis(cfg, comp);

% none of these should fail
ft_checkdata(freq, 'datatype', 'freq');
ft_checkdata(freq, 'datatype', 'comp');
ft_checkdata(freq, 'datatype', 'freq+comp');
