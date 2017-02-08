function test_bug3218

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug3218
% TEST ft_rejectartifact
% TEST ft_rejectvisual

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3218.mat'), 'data');

% When rejecting artifacts using 'partial', ft_rejectvisual works fine
cfg = [];
cfg.artfctdef.visual.artifact = data.sampleinfo(1, :) + [300 -300];
cfg.artfctdef.reject = 'partial';
dataReject = ft_rejectartifact(cfg, data);

ft_rejectvisual([], dataReject);

% When rejecting artifacts using 'nan', ft_rejectvisual recalculates the
% metric after each trial that is excluded. With this small data set it's
% not really a problem, but it becomes cumbersome with many trials and
% channels.

cfg = [];
cfg.artfctdef.visual.artifact = data.sampleinfo(1, :) + [300 -300];
cfg.artfctdef.reject = 'nan';
dataReject = ft_rejectartifact(cfg, data);

ft_rejectvisual([], dataReject);


