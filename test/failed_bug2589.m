function failed_bug2589

% WALLTIME 00:10:00
% MEM 1gb

% TEST test_bug2589
% TEST ft_selectdata getdimord

data = [];
for i=1:2
  data.label{i} = num2str(i);
end

for i=1:7
  data.trial{i} = randn(2,300);
  data.time{i}  = (1:300)/300;
end

data.sampleinfo = randn(7, 2); % 7 trials, not 2 channels
data.trialinfo  = randn(7, 7); % 7 trials, 7 "triggers"

cfg = [];
cfg.channel = 1;
output = ft_selectdata(cfg, data);

assert(size(output.sampleinfo,1)==7);
assert(size(output.sampleinfo,2)==2); % this failed

%% do a similar check on trial selections
cfg = [];
cfg.trials = [2 3 4];
output = ft_selectdata(cfg, data);

assert(size(output.sampleinfo,1)==3);
assert(size(output.sampleinfo,2)==2);
assert(size(output.trialinfo,1)==3);
assert(size(output.trialinfo,2)==7);


