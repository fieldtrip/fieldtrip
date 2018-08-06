function test_bug27

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preproc_lowpassfilter ft_resampledata ft_resampledata 

% this script tests for bug 27 and for the solution
% the 'bug' is that if ft_resampledata is called on data with different time axes
% after resampling the time axes of the trials may be shifted with respect to one
% another with a fractional amount (governed by the input fsample)

% generate some data
data          = [];
data.trial{1} = randn(2,1000);
data.time{1}  = [0:0.001:0.999] - 0.5;
data.trial{2} = randn(2,1000);
data.time{2}  = [0.001:0.001:1] - 0.5;
data.label    = {'chan1';'chan2'};
data.fsample  = 1000;
data.cfg.trl  = zeros(2,3);

for k = 1:2
  data.trial{k} = ft_preproc_lowpassfilter(data.trial{k}, 1000, 40);
end

% resample
cfg            = [];
cfg.resamplefs = 256;
cfg.detrend    = 'no';
data2          = ft_resampledata(cfg, data);

% make it a bit more extreme
data.time{2} = data.time{2} + 3.1152;
data3        = ft_resampledata(cfg, data);

