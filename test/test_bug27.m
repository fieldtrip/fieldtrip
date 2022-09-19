function test_bug27

% MEM 8gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preproc_lowpassfilter ft_resampledata ft_resampledata ft_datatype_raw fixtimeaxes

% this script tests for bug 27 and for the solution
% the 'bug' is that if ft_resampledata is called on data with different time axes
% after resampling the time axes of the trials may be shifted with respect to one
% another with a fractional amount (governed by the input fsample)

%% generate some data
% shifted by one sample
data1          = [];
data1.trial{1} = randn(2,1000);
data1.trial{2} = randn(2,1000);
data1.time{1}  = [0:0.001:0.999] - 0.5;
data1.time{2}  = [0.001:0.001:1] - 0.5;
data1.label    = {'chan1';'chan2'};
data1.fsample  = 1000;
data1.cfg.trl  = zeros(2,3);

for k = 1:2
  data1.trial{k} = ft_preproc_lowpassfilter(data1.trial{k}, 1000, 40);
end


%% resample with a factor 1x
cfg            = [];
cfg.resamplefs = 1000;
cfg.detrend    = 'no';
data2          = ft_resampledata(cfg, data1);

[ok, msg] = isalmostequal(data1.time, data2.time, 'abstol', 1e-10);
assert(ok, msg);


%% resample
cfg            = [];
cfg.resamplefs = 256;
cfg.detrend    = 'no';
data2          = ft_resampledata(cfg, data1);

figure
subplot(2,2,1), plot(data1.time{1}, data1.trial{1}, 'b'); xlim([data1.time{1}(1) data1.time{1}(end)])
subplot(2,2,2), plot(data1.time{2}, data1.trial{2}, 'b'); xlim([data1.time{2}(1) data1.time{2}(end)])
subplot(2,2,3), plot(data2.time{1}, data2.trial{1}, 'r'); xlim([data1.time{1}(1) data1.time{1}(end)])
subplot(2,2,4), plot(data2.time{2}, data2.trial{2}, 'r'); xlim([data1.time{2}(1) data1.time{2}(end)])


%% make it a bit more extreme
data1.time{2} = data1.time{2} + 3.1152;
data2 = ft_resampledata(cfg, data1);

figure
subplot(2,2,1), plot(data1.time{1}, data1.trial{1}, 'b'); xlim([data1.time{1}(1) data1.time{1}(end)])
subplot(2,2,2), plot(data1.time{2}, data1.trial{2}, 'b'); xlim([data1.time{2}(1) data1.time{2}(end)])
subplot(2,2,3), plot(data2.time{1}, data2.trial{1}, 'r'); xlim([data1.time{1}(1) data1.time{1}(end)])
subplot(2,2,4), plot(data2.time{2}, data2.trial{2}, 'r'); xlim([data1.time{2}(1) data1.time{2}(end)])


%% make it a bit more extreme
cfg.resamplefs = 17;
data2 = ft_resampledata(cfg, data1);

figure
subplot(2,2,1), plot(data1.time{1}, data1.trial{1}, 'b'); xlim([data1.time{1}(1) data1.time{1}(end)])
subplot(2,2,2), plot(data1.time{2}, data1.trial{2}, 'b'); xlim([data1.time{2}(1) data1.time{2}(end)])
subplot(2,2,3), plot(data2.time{1}, data2.trial{1}, 'r'); xlim([data1.time{1}(1) data1.time{1}(end)])
subplot(2,2,4), plot(data2.time{2}, data2.trial{2}, 'r'); xlim([data1.time{2}(1) data1.time{2}(end)])


%%

% this contains a problematic data structure from Mats
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug27.mat'), 'data');

clear sample0 time0
for i=1:numel(data.time)
  % find the sample corresponding with time zero
  sample0(i)  = nearest(data.time{i}, 0);
  time0(i)    = data.time{i}(sample0(i));
end
figure
plot(time0, '.-'); % this looks ok, no jitter
assert(std(time0)<1e-6, 'there is too much jitter in time0');

cfg = [];
cfg.resamplefs = 600;
data_resampled = ft_resampledata(cfg, data);

clear sample0 time0
for i=1:numel(data.time)
  % find the sample corresponding with time zero
  sample0(i)  = nearest(data_resampled.time{i}, 0);
  time0(i)    = data_resampled.time{i}(sample0(i));
end
figure
plot(time0, '.-'); % this showed quite some jitter over trials
assert(std(time0)<1e-6, 'there is too much jitter in time0');

