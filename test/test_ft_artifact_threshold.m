function test_ft_artifact_threshold

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY

% this script serves to test the asymetric onset/offset and the peak-detection

%%

data = [];
data.label = {'ramp'};
data.fsample = 1000;
data.time{1} = (1:1000)/1000;
data.trial{1} = zeros(1,1000);
data.sampleinfo = [1 1000];

% ramp up and down
data.trial{1}(  1: 500) = 1:500;
data.trial{1}(501:1000) = 500:-1:1;

%%

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.max = 400;
[cfg1, artifact1] = ft_artifact_threshold(cfg, data);

cfg = [];
cfg.trl = artifact1;
data1 = ft_redefinetrial(cfg, data);

cfg = [];
timelock1 = ft_timelockanalysis(cfg, data1);

figure
plot(timelock1.time, timelock1.avg, '.-');

% the peak should be at t=0
[val, indx] = max(timelock1.avg);
assert(isalmostequal(timelock1.time(indx), 0));

%%

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.onset  = 400;
cfg.artfctdef.threshold.offset = 400;
[cfg2, artifact2] = ft_artifact_threshold(cfg, data);

assert(isequal(artifact2, artifact1));

%%

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.onset  = 400;
cfg.artfctdef.threshold.offset = 300;
[cfg3, artifact3] = ft_artifact_threshold(cfg, data);

assert(~isequal(artifact3, artifact1));

cfg = [];
cfg.trl = artifact3;
data3 = ft_redefinetrial(cfg, data);

cfg = [];
timelock3 = ft_timelockanalysis(cfg, data3);

figure
plot(timelock3.time, timelock3.avg, '.-'); % it should be asymetric

% the peak should be at t=0
[val, indx] = max(timelock3.avg);
assert(isalmostequal(timelock3.time(indx), 0, 'abstol', 100*eps));


%%

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.onset  = 1;
cfg.artfctdef.threshold.offset = 1;
[cfg4, artifact4] = ft_artifact_threshold(cfg, data);

assert(isequal(artifact4, [1 1000 -499]));

%%
% invert the data and invert the thresholds

data.trial{1} = -data.trial{1};

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.onset  = -400;
cfg.artfctdef.threshold.offset = -300;
[cfg5, artifact5] = ft_artifact_threshold(cfg, data);

assert(isequal(artifact5, artifact3));

