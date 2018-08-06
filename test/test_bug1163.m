function test_bug1163

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1163 ft_timelockbaseline

% the problem seems to be that ft_timelockbaseline explicitly only checks
% for a avg field

% Setup
data.label = {'a', 'b', 'c', 'd'};
data.time = linspace(-1, 1, 400);
data.individual = zeros(19, 4, 400);
data.dimord = 'subj_chan_time';
data.cfg = {};

% Example given in bug 1163
cfg = [];
cfg.baseline = [-0.1 0];
data = ft_timelockbaseline(cfg, data) % produces error;
