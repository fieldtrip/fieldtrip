function test_bug1245

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug_1245
% TEST ft_multiplotER

% The issue: when inputting data where the corresponding layout consists of
% more channels than the to-be-plotted channels (specified in cfg.channel),
% and where those channels may have a different order of magnitude, the
% automatic scaling of the y-axis is incorrect, because the scale
% determination algorithm includes the channels that are excluded from the
% plotting

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))

% this is one way of having ft_multiplotER select avg as the parameter
load bug1245.mat
data.avg = [];
figure
ft_multiplotER(cfg, data);

% this is another one way of having ft_multiplotER select avg as the parameter
load bug1245.mat
cfg.parameter = 'avg';
figure
ft_multiplotER(cfg, data);

% but the data contains a 'trial' field, so better explicitly select that
load bug1245.mat
cfg.parameter = 'trial';
figure
ft_multiplotER(cfg, data);
