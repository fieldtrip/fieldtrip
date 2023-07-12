function test_issue2072

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_selectdata
% DATA no
 
% Simulate data
data = [];
data.fsample = 1000;
data.label   = {'1', '2', '3'};
for i = 1:3
    data.time{i}  = -0.5:1/data.fsample:i;
    data.trial{i} = randn(length(data.label), length(data.time{i}));
end

% This originally cut the TOI to [-0.5, 1], which is unexpected, given that
% only the second trial is selected, which has a long enough time axis
cfg = [];
cfg.trials  = 2;
cfg.latency = [-0.5, 1.5];
data_cut    = ft_selectdata(cfg, data);
assert(data_cut.time{1}(end)==1.5);

% This will give the correct TOI of [-0.5, 1.5]
cfg = [];
cfg.trials  = 2;
data_cut    = ft_selectdata(cfg, data);

cfg = [];
cfg.latency = [-0.5, 1.5];
data_cut    = ft_selectdata(cfg, data_cut);
assert(data_cut.time{1}(end)==1.5);
