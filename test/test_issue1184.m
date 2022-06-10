function test_issue1184

% MEM 8gb
% WALLTIME 00:30:00
% DEPENDENCY ft_resampledata

data = [];
data.time{1} = (1:228864) * 1/1200;
data.trial{1}(1,:) = 1:228864;
data.label = {'time'};
data.fsample = 1200;

cfg = [];
cfg.resamplefs = 100;
resampled = ft_resampledata(cfg, data);

assert(length(resampled.time{1})==19072);

%%

data11 = [];
data11.time{1} = (1:11) * 1/1200;
data11.trial{1}(1,:) = 1:11;
data11.label = {'time'};
data11.fsample = 1200;

cfg = [];
cfg.resamplefs = 100;
resampled11 = ft_resampledata(cfg, data11);

assert(isalmostequal(resampled11.time{1}(1), mean(data11.time{1}), 'abstol', 1e-6));

data12 = [];
data12.time{1} = (1:12) * 1/1200;
data12.trial{1}(1,:) = 1:12;
data12.label = {'time'};
data12.fsample = 1200;

cfg = [];
cfg.resamplefs = 100;
resampled12 = ft_resampledata(cfg, data12);

assert(isalmostequal(resampled12.time{1}(1), mean(data12.time{1}), 'abstol', 1e-6));


%%

% resampled11 and resampled12 have a time axis that is shifted by half a sample (at 1200 Hz)
% since they are constructed from two data structures that are different by one sample
%
% the following code would result in a data structure in which the trials are slightly mis-aligned
% resampledXX = ft_appenddata([], resampled11, resampled12);

% when resampling both trials at the same time, the time axes should remain consistent
dataXX = ft_appenddata([], data11, data12);

for resamplefs=1:1200
  cfg = [];
  cfg.showcallinfo = 'no';
  cfg.resamplefs = resamplefs;
  resampledXX = ft_resampledata(cfg, dataXX);
  
  % with a resampling rate of 600 or higher, the time axes of the resampled trials will have a different length due to the padding
  % only the first part should be compared
  n = min(length(resampledXX.time{1}), length(resampledXX.time{2}));
  assert(isalmostequal(mean(resampledXX.time{1}(1:n)), mean(resampledXX.time{2}(1:n)), 'abstol', 1e-6));
end

