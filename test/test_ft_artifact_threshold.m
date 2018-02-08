function test_ft_artifact_threshold

% WALLTIME 00:10:00
% MEM 1gb

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

%%

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.onset  = 1;
cfg.artfctdef.threshold.offset = 1;
[cfg4, artifact4] = ft_artifact_threshold(cfg, data);

assert(isequal(artifact4, [1 1000]));

%%
% invert the data and invert the thresholds

data.trial{1} = -data.trial{1};

cfg = [];
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.onset  = -400;
cfg.artfctdef.threshold.offset = -300;
[cfg5, artifact5] = ft_artifact_threshold(cfg, data);

assert(isequal(artifact5, artifact3));


